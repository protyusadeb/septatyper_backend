"""Background analysis, flagging and curator approval for community submissions.

Nothing here touches the live tables (isolates, septatypes, ...) until a curator approves.
"""
from __future__ import annotations

import gzip
import logging
import os
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import requests
from sqlalchemy import func, select, update
from sqlalchemy.orm import Session

from .bvbrc import BVBRCClient, summarise_computational
from .engine_runtime import EngineBusy, run_engine
from .models import (
    AuditLog,
    BvbrcPrediction,
    Habitat,
    Isolate,
    IsolateGenotype,
    IsolateSource,
    Phenotype,
    PMQRCall,
    SeptaPrediction,
    Septatype,
    Species,
    Submission,
    SubmissionStatus,
    User,
)
from .prediction_engine import ECOLI, SSONNEI, EngineError, InputError, aa_level_septatype

log = logging.getLogger("septatyper.submissions")

BVBRC_API_BASE = os.getenv("BVBRC_API_BASE", "https://www.bv-brc.org/api")

# code -> (label, severity). "review" flags put the submission in the needs_review state
# and require the curator to confirm before approving. "info" flags are shown only.
FLAGS: Dict[str, tuple] = {
    "novel_septatype": ("Novel septatype (not in the database)", "review"),
    "phenotype_mismatch": ("Predicted phenotype conflicts with the submitted phenotype", "review"),
    "indeterminate_prediction": ("Prediction could not be completed", "review"),
    "species_mismatch": ("Declared species differs from the species inferred from the genome", "review"),
    "bvbrc_id_not_found": ("BV-BRC genome ID was not found", "review"),
    "bvbrc_lookup_failed": ("BV-BRC lookup failed (service error)", "review"),
    "duplicate_genome": ("The same genome sequence was already submitted", "review"),
    "mic_phenotype_conflict": ("MIC contradicts the submitted phenotype", "review"),
    "new_codon_variant": ("Known septatype with a new RecN codon notation", "info"),
}
REVIEW_FLAGS = {code for code, (_, sev) in FLAGS.items() if sev == "review"}


def flag_objects(codes: Optional[List[str]]) -> List[dict]:
    return [{"code": c, "label": FLAGS[c][0], "severity": FLAGS[c][1]} for c in (codes or []) if c in FLAGS]


def species_group(name: Optional[str]) -> Optional[str]:
    """E. coli and S. sonnei cannot be separated by gyrA/parC/recN, so they count as one group."""
    if name in (ECOLI, SSONNEI):
        return "Escherichia/Shigella"
    return name


def mic_conflicts_with_phenotype(mic: Optional[float], unit: Optional[str], phenotype: str) -> bool:
    """True only where CLSI and EUCAST agree: MIC >= 1 mg/L is never Susceptible, MIC <= 0.25 never Resistant."""
    if mic is None or unit not in ("mg/L", "ug/mL"):
        return False
    return (mic >= 1.0 and phenotype == "Susceptible") or (mic <= 0.25 and phenotype == "Resistant")


def _now() -> datetime:
    return datetime.now(timezone.utc)


# ===================================================================== background work
def process_submission(bind, submission_id: int) -> None:
    """Runs after the upload response is sent: BV-BRC lookup, prediction engine, flags."""
    with Session(bind) as db:
        sub = db.get(Submission, submission_id)
        if sub is None:
            return
        try:
            sub.analysis_state = "running"
            db.commit()
            _lookup_bvbrc(sub)
            result = _run_prediction(sub)
            _evaluate(db, sub, result)
            sub.analysis_state = "done"
            sub.analysis_error = None
        except InputError as exc:
            sub.analysis_state, sub.analysis_error = "failed", f"The genome could not be analysed: {exc}"
        except EngineBusy:
            sub.analysis_state, sub.analysis_error = "failed", "The analysis server was busy. Please resubmit."
        except EngineError:
            log.exception("Engine failure for submission %s", submission_id)
            sub.analysis_state, sub.analysis_error = "failed", "The analysis pipeline failed. A curator will look into it."
        except Exception:  # never leave a submission stuck in 'running'
            log.exception("Unexpected failure for submission %s", submission_id)
            db.rollback()
            sub = db.get(Submission, submission_id)
            sub.analysis_state, sub.analysis_error = "failed", "Unexpected error during analysis."
        db.commit()


def _lookup_bvbrc(sub: Submission) -> None:
    if not sub.bv_brc_id_supplied:
        sub.bv_brc_lookup_status = "not_provided"
        return
    client = BVBRCClient(base=BVBRC_API_BASE, sleep=0.2, batch_size=1, timeout=15, retries=2)
    try:
        rec = client.fetch_batch([sub.bv_brc_id_supplied])[sub.bv_brc_id_supplied]
    except (requests.RequestException, ValueError, KeyError) as exc:
        log.warning("BV-BRC lookup failed for %s: %s", sub.bv_brc_id_supplied, exc)
        sub.bv_brc_lookup_status = "error"
        return
    name = rec.get("genome_name")
    if not name:
        sub.bv_brc_lookup_status = "not_found"
        return
    pred, raw, method = summarise_computational(rec.get("amr", []))
    sub.bv_brc_name_fetched = name
    sub.bv_brc_lookup_status = "ok"
    sub.bvbrc_amr = {"prediction": pred.value, "raw": raw, "method": method}


def _run_prediction(sub: Submission) -> dict:
    with gzip.open(sub.file_path, "rt", encoding="utf-8") as fh:  # type: ignore[arg-type]
        text = fh.read()
    return run_engine(text, sub.declared_species, wait=300.0)


def _evaluate(db: Session, sub: Submission, result: dict) -> None:
    st: Optional[str] = result["septatype"]
    pred: str = result["prediction"]
    flags: List[str] = []

    sub.prediction_result = result
    sub.predicted_septatype_string = st
    sub.septatyper_prediction = pred

    exact = db.scalar(select(Septatype).where(Septatype.septatype_string == st)) if st else None
    plain = aa_level_septatype(st)
    known = exact or (db.scalar(select(Septatype).where(Septatype.septatype_string == plain)) if plain else None)
    sub.predicted_septatype_id = exact.septatype_id if exact else None
    sub.is_novel_septatype = bool(st) and known is None

    if pred == "Indeterminate" or st is None:
        flags.append("indeterminate_prediction")
        sub.prediction_matches_phenotype = None
    else:
        if known is None:
            flags.append("novel_septatype")
        elif exact is None:
            flags.append("new_codon_variant")
        expected = "Cres" if sub.experimental_phenotype == Phenotype.RESISTANT else "Csus"
        sub.prediction_matches_phenotype = pred == expected
        if not sub.prediction_matches_phenotype:
            flags.append("phenotype_mismatch")

    inferred = result["species"].get("inferred")
    if inferred and species_group(inferred) != species_group(sub.declared_species):
        flags.append("species_mismatch")
    if sub.bv_brc_lookup_status == "not_found":
        flags.append("bvbrc_id_not_found")
    elif sub.bv_brc_lookup_status == "error":
        flags.append("bvbrc_lookup_failed")
    if sub.file_sha256 and db.scalar(
        select(func.count()).select_from(Submission).where(
            Submission.file_sha256 == sub.file_sha256,
            Submission.submission_id != sub.submission_id,
            Submission.status != SubmissionStatus.REJECTED,
        )
    ):
        flags.append("duplicate_genome")
    if sub.experimental_phenotype and mic_conflicts_with_phenotype(
        float(sub.experimental_mic) if sub.experimental_mic is not None else None,
        sub.mic_unit, sub.experimental_phenotype.value,
    ):
        flags.append("mic_phenotype_conflict")

    sub.flags = flags
    sub.status = (
        SubmissionStatus.NEEDS_REVIEW if REVIEW_FLAGS.intersection(flags) else SubmissionStatus.PENDING
    )


# ======================================================================= curator actions
class ApprovalError(Exception):
    def __init__(self, status_code: int, message: str):
        super().__init__(message)
        self.status_code, self.message = status_code, message


def review_flags(sub: Submission) -> List[str]:
    return [f for f in (sub.flags or []) if f in REVIEW_FLAGS]


def _check_reviewable(sub: Optional[Submission]) -> Submission:
    if sub is None:
        raise ApprovalError(404, "Submission not found")
    if sub.status not in (SubmissionStatus.PENDING, SubmissionStatus.NEEDS_REVIEW):
        raise ApprovalError(409, f"Submission is already {sub.status.value}")
    if sub.analysis_state not in ("done", "failed"):
        raise ApprovalError(409, "Analysis is still running")
    return sub


def approve_submission(db: Session, submission_id: int, curator: User, notes: Optional[str], confirm_flags: bool) -> dict:
    """Copy a reviewed submission into the live tables. The caller commits (one transaction)."""
    sub = _check_reviewable(db.get(Submission, submission_id, with_for_update=True))
    if sub.analysis_state != "done" or not sub.prediction_result or not sub.predicted_septatype_string:
        raise ApprovalError(409, "This submission has no completed prediction and can only be rejected")
    flags = review_flags(sub)
    if flags and not confirm_flags:
        raise ApprovalError(409, "Flagged submission: confirm the flags to approve (confirm_flags=true). Flags: " + ", ".join(flags))

    result = sub.prediction_result
    st = sub.predicted_septatype_string
    genome_id = (
        sub.bv_brc_id_supplied if sub.bv_brc_lookup_status == "ok" and sub.bv_brc_id_supplied
        else f"SEPTA-{sub.submission_id:06d}"
    )
    if db.get(Isolate, genome_id) is not None:
        raise ApprovalError(409, f"Genome ID {genome_id} is already in the database")

    species = db.scalar(select(Species).where(Species.scientific_name == sub.declared_species))
    if species is None:
        species = Species(scientific_name=sub.declared_species)
        db.add(species)
        db.flush()
    habitat = None
    if sub.habitat:
        habitat = db.scalar(select(Habitat).where(Habitat.name == sub.habitat))
        if habitat is None:
            habitat = Habitat(name=sub.habitat)
            db.add(habitat)
            db.flush()

    codons = result["codons"]
    aa = {k: c["aa"] for k, c in codons.items()}
    septa = db.scalar(select(Septatype).where(Septatype.septatype_string == st))
    created_septatype = septa is None
    if septa is None:
        septa = Septatype(
            septatype_string=st,
            gyra83_aa=aa["gyrA_83"], gyra87_aa=aa["gyrA_87"], parc80_aa=aa["parC_80"],
            parc84_aa=aa["parC_84"], recn346_aa=aa["recN_346"],
            qnrb2=result["pmqr"]["qnrB2"]["present"], qnrs1=result["pmqr"]["qnrS1"]["present"],
            mean_hydrophobicity=result.get("mean_hydrophobicity"),
            is_novel=bool(sub.is_novel_septatype), first_seen_date=date.today(), isolate_count=0,
        )
        db.add(septa)
        db.flush()

    bv = sub.bvbrc_amr if sub.bv_brc_lookup_status == "ok" and sub.bvbrc_amr else None
    pred = sub.septatyper_prediction
    year = sub.collection_year
    isolate = Isolate(
        bv_brc_genome_id=genome_id,
        official_genome_name=sub.bv_brc_name_fetched or sub.isolate_name or f"Community submission {sub.submission_id}",
        isolate_name=sub.isolate_name,
        species_id=species.species_id,
        habitat_id=habitat.habitat_id if habitat else None,
        collection_year_start=year, collection_year_end=year,
        collection_year_raw=str(year) if year else None,
        country_name=sub.country, country_iso3=sub.country_iso3,
        lab_phenotype=sub.experimental_phenotype,
        mic_value=sub.experimental_mic, mic_unit=sub.mic_unit, mic_method=sub.mic_method,
        bv_brc_computational_prediction=BvbrcPrediction(bv["prediction"]) if bv else None,
        bv_brc_prediction_raw=bv["raw"] if bv else None,
        bv_brc_prediction_method=bv["method"] if bv else None,
        bv_brc_fetched_at=_now() if bv else None,
        septatyper_prediction=SeptaPrediction(pred) if pred in ("Cres", "Csus") else None,
        septatype_id=septa.septatype_id,
        source=IsolateSource.COMMUNITY,
        submission_id=sub.submission_id,
    )
    db.add(isolate)
    db.flush()

    genes = [g for g in result["genes"].values() if g.get("found")]
    db.add(IsolateGenotype(
        bv_brc_genome_id=genome_id,
        gyra83_codon=codons["gyrA_83"]["codon"], gyra83_aa=aa["gyrA_83"],
        gyra87_codon=codons["gyrA_87"]["codon"], gyra87_aa=aa["gyrA_87"],
        parc80_codon=codons["parC_80"]["codon"], parc80_aa=aa["parC_80"],
        parc84_codon=codons["parC_84"]["codon"], parc84_aa=aa["parC_84"],
        recn346_codon=codons["recN_346"]["codon"], recn346_aa=aa["recN_346"],
        recn346_synonymous_flag=result["recn346"]["synonymous_flag"],
        alignment_identity=min((g["identity"] for g in genes), default=None),
        coverage=min((g["coverage"] for g in genes), default=None),
        reference_used=";".join(f"{k}:{v['reference']}" for k, v in result["genes"].items() if v.get("found"))[:80],
        pipeline_version=result["pipeline_version"],
    ))
    for gene, call in result["pmqr"].items():
        db.add(PMQRCall(
            bv_brc_genome_id=genome_id, gene=gene, present=bool(call["present"]),
            identity=call["identity"], coverage=call["coverage"], contig_id=call["contig"],
        ))
    db.flush()
    db.execute(
        update(Septatype).where(Septatype.septatype_id == septa.septatype_id).values(
            isolate_count=select(func.count(Isolate.bv_brc_genome_id))
            .where(Isolate.septatype_id == septa.septatype_id).scalar_subquery()
        )
    )

    sub.status = SubmissionStatus.APPROVED
    sub.curator_id, sub.curator_notes, sub.decided_at = curator.user_id, notes, _now()
    db.add(AuditLog(
        user_id=curator.user_id, action="submission_approved", target_table="submissions",
        target_id=str(sub.submission_id),
        detail={"isolate": genome_id, "septatype": st, "created_septatype": created_septatype,
                "confirmed_flags": flags, "notes": notes},
    ))
    return {"isolate_id": genome_id, "septatype_id": septa.septatype_id,
            "septatype_string": st, "created_septatype": created_septatype}


def reject_submission(db: Session, submission_id: int, curator: User, notes: str) -> None:
    sub = _check_reviewable(db.get(Submission, submission_id, with_for_update=True))
    sub.status = SubmissionStatus.REJECTED
    sub.curator_id, sub.curator_notes, sub.decided_at = curator.user_id, notes, _now()
    db.add(AuditLog(
        user_id=curator.user_id, action="submission_rejected", target_table="submissions",
        target_id=str(sub.submission_id), detail={"notes": notes},
    ))
