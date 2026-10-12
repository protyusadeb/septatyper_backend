"""Ingest Table S1 into the SeptaTyper database.

Usage
    python -m septatyper.ingest --xlsx "Table S1.xlsx"
    python -m septatyper.ingest --xlsx "Table S1.xlsx" --dry-run
    python -m septatyper.ingest --xlsx "Table S1.xlsx" --skip-api --limit 50

What it does
1. Parses the sheet. Species and habitat come from the 7 section header rows
   (e.g. "Escherichia coli FECAL"), not from columns.
2. Cleans years (single, range, invalid), MLST nulls, countries, qnr flags.
3. Queries BV-BRC for the official genome_name and the ciprofloxacin
   computational phenotype (genome_amr collection, evidence = Computational Method).
   Responses are cached to a JSON file so an interrupted run resumes cheaply.
4. Upserts species, habitats, septatypes, isolates, genotypes (amino acids only)
   and pmqr_calls. Safe to re-run: existing codon data is never overwritten.
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
import time
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from sqlalchemy import func, select, update
from sqlalchemy.orm import Session
from urllib3.util.retry import Retry

from .bvbrc import BVBRCClient, summarise_computational
from .countries import country_to_iso3
from .database import SessionLocal, get_engine, init_db
from .models import (
    BvbrcPrediction,
    Habitat,
    Isolate,
    IsolateGenotype,
    IsolateSource,
    Phenotype,
    PMQRCall,
    Septatype,
    Species,
)
from .seed import AMINO_ACIDS, seed_lookups

log = logging.getLogger("septatyper.ingest")

EXPECTED_COLUMNS = [
    "Genome ID", "Isolate Name", "Ciprofloxacin Resistant Phenotype", "Collection Year",
    "Isolation Country", "MLST", "GyrA-83", "GyrA-87", "ParC-80", "ParC-84",
    "RecN-346", "qnrB2", "qnrS1",
]
GENOME_ID_RE = re.compile(r"^\d+\.\d+$")
SECTION_RE = re.compile(r"^(?P<species>.+?)\s+(?P<habitat>[A-Za-z]+)$")
YEAR_RE = re.compile(r"^(\d{4})(?:\s*[-\u2013]\s*(\d{4}))?$")
AA_SET = set(AMINO_ACIDS)


# =============================================================================
# 1. Parsing and cleaning
# =============================================================================
def _clean(value: Any) -> Optional[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    s = str(value).strip()
    return None if s == "" or s.lower() in {"nan", "none", "na", "n/a"} else s


def parse_year(raw: Optional[str], warnings: List[str], genome_id: str) -> Tuple[Optional[int], Optional[int], Optional[str]]:
    """Return (start, end, raw). Single year -> start == end. Invalid parts -> None."""
    s = _clean(raw)
    if s is None:
        return None, None, None
    s = re.sub(r"\.0$", "", s)
    m = YEAR_RE.match(s)
    if not m:
        warnings.append(f"{genome_id}: unparseable year '{s}'")
        return None, None, s
    first, second = int(m.group(1)), (int(m.group(2)) if m.group(2) else None)
    this_year = date.today().year

    def valid(y: Optional[int]) -> Optional[int]:
        return y if y is not None and 1900 <= y <= this_year else None

    start, end = valid(first), valid(second if second is not None else first)
    if (start is None) or (second is not None and end is None):
        warnings.append(f"{genome_id}: implausible year '{s}' (start={start}, end={end}); raw value kept")
    return start, end, s


def _isna(v: Any) -> bool:
    return v is None or (isinstance(v, float) and pd.isna(v))


def parse_binary(value: Any) -> Optional[int]:
    s = _clean(value)
    if s is None:
        return None
    try:
        v = int(float(s))
    except ValueError:
        return None
    return v if v in (0, 1) else None


def parse_table_s1(path: Path) -> Tuple[pd.DataFrame, List[str]]:
    """Return (tidy dataframe, warnings)."""
    warnings: List[str] = []
    raw = pd.read_excel(path, header=None, dtype=str)

    hdr_rows = raw.index[raw[0].astype(str).str.strip().eq("Genome ID")]
    if len(hdr_rows) == 0:
        raise ValueError("Header row 'Genome ID' not found in column A")
    hdr_i = hdr_rows[0]
    header = [_clean(x) for x in raw.loc[hdr_i].tolist()][: len(EXPECTED_COLUMNS)]
    if header != EXPECTED_COLUMNS:
        raise ValueError(f"Unexpected header.\n expected: {EXPECTED_COLUMNS}\n found:    {header}")

    species: Optional[str] = None
    habitat: Optional[str] = None
    records: List[Dict[str, Any]] = []
    sections: List[str] = []

    for _, row in raw.loc[hdr_i + 1:].iterrows():
        cells = [_clean(x) for x in row.tolist()[: len(EXPECTED_COLUMNS)]]
        if cells[0] is None:
            continue
        if all(c is None for c in cells[1:]):  # section header row
            m = SECTION_RE.match(cells[0])
            if not m:
                raise ValueError(f"Cannot split section header into species + habitat: '{cells[0]}'")
            species, habitat = m.group("species"), m.group("habitat").capitalize()
            sections.append(cells[0])
            continue
        if species is None:
            raise ValueError(f"Isolate row before any section header: {cells[0]}")
        rec = dict(zip(EXPECTED_COLUMNS, cells))
        rec["species"], rec["habitat"] = species, habitat
        records.append(rec)

    df = pd.DataFrame(records)
    log.info("Found %d section headers: %s", len(sections), sections)

    # Genome ID checks
    bad_ids = df[~df["Genome ID"].str.match(GENOME_ID_RE)]
    if len(bad_ids):
        raise ValueError(f"{len(bad_ids)} malformed Genome IDs, e.g. {bad_ids['Genome ID'].head().tolist()}")
    dupes = df["Genome ID"][df["Genome ID"].duplicated()].unique().tolist()
    if dupes:
        raise ValueError(f"Duplicate Genome IDs: {dupes[:10]}")

    # Phenotype
    pheno_map = {"resistant": Phenotype.RESISTANT, "susceptible": Phenotype.SUSCEPTIBLE}
    df["lab_phenotype"] = df["Ciprofloxacin Resistant Phenotype"].str.lower().map(pheno_map)
    bad = df[df["lab_phenotype"].isna()]
    if len(bad):
        raise ValueError(f"Unrecognised phenotype values: {bad['Ciprofloxacin Resistant Phenotype'].unique().tolist()}")

    # Years
    years = [parse_year(r, warnings, g) for r, g in zip(df["Collection Year"], df["Genome ID"])]
    df["year_start"], df["year_end"], df["year_raw"] = zip(*years) if years else ([], [], [])

    # Country, MLST
    df["country_iso3"] = df["Isolation Country"].map(country_to_iso3)
    for g, c in df.loc[df["country_iso3"].isna() & df["Isolation Country"].notna(), ["Genome ID", "Isolation Country"]].itertuples(index=False):
        warnings.append(f"{g}: no ISO3 code for country '{c}'")
    df["MLST"] = df["MLST"].map(lambda v: None if _isna(v) else _clean(v))
    df["MLST"] = df["MLST"].map(lambda v: None if _isna(v) else (v if v.upper().startswith("ST") else f"ST{v}"))

    # Amino acids
    for col in ["GyrA-83", "GyrA-87", "ParC-80", "ParC-84", "RecN-346"]:
        df[col] = df[col].map(lambda v: None if _isna(v) else v.upper())
        invalid = df[df[col].notna() & ~df[col].isin(AA_SET)]
        for g, v in invalid[["Genome ID", col]].itertuples(index=False):
            warnings.append(f"{g}: invalid amino acid '{v}' in {col}; set to NULL")
        df.loc[~df[col].isin(AA_SET), col] = None

    # PMQR flags
    for col in ["qnrB2", "qnrS1"]:
        df[col] = df[col].map(parse_binary)
        n_bad = int(df[col].isna().sum())
        if n_bad:
            warnings.append(f"{n_bad} rows with missing/invalid {col}")

    # Septatype string
    def septa(r: pd.Series) -> Optional[str]:
        parts = [r["GyrA-83"], r["GyrA-87"], r["ParC-80"], r["ParC-84"], r["RecN-346"], r["qnrB2"], r["qnrS1"]]
        return None if any(p is None or (isinstance(p, float) and pd.isna(p)) for p in parts) else "-".join(str(int(p)) if isinstance(p, (int, float)) else p for p in parts)

    df["septatype_string"] = df.apply(septa, axis=1)
    df = df.astype(object).where(df.notna(), None)
    return df, warnings


# BV-BRC client and computational-phenotype summary live in septatyper/bvbrc.py


# =============================================================================
# 3. Database load
# =============================================================================
def _get_or_create_species(session: Session, cache: Dict[str, int], name: str) -> int:
    if name not in cache:
        sp = session.scalar(select(Species).where(Species.scientific_name == name))
        if sp is None:
            sp = Species(scientific_name=name)
            session.add(sp)
            session.flush()
            log.warning("Species '%s' was not in the seed list; created without taxid", name)
        cache[name] = sp.species_id
    return cache[name]


def _get_or_create_habitat(session: Session, cache: Dict[str, int], name: Optional[str]) -> Optional[int]:
    if name is None:
        return None
    if name not in cache:
        h = session.scalar(select(Habitat).where(Habitat.name == name))
        if h is None:
            h = Habitat(name=name)
            session.add(h)
            session.flush()
        cache[name] = h.habitat_id
    return cache[name]


def _get_or_create_septatype(session: Session, cache: Dict[str, int], row: Dict[str, Any]) -> Optional[int]:
    s = row["septatype_string"]
    if s is None:
        return None
    if s not in cache:
        st = Septatype(
            septatype_string=s,
            gyra83_aa=row["GyrA-83"], gyra87_aa=row["GyrA-87"], parc80_aa=row["ParC-80"],
            parc84_aa=row["ParC-84"], recn346_aa=row["RecN-346"],
            qnrb2=int(row["qnrB2"]), qnrs1=int(row["qnrS1"]),
            first_seen_date=date.today(),
        )
        session.add(st)
        session.flush()
        cache[s] = st.septatype_id
    return cache[s]


def load_isolates(session: Session, df: pd.DataFrame, bvbrc: Optional[Dict[str, Dict[str, Any]]]) -> Counter:
    stats: Counter = Counter()
    sp_cache: Dict[str, int] = {}
    hab_cache: Dict[str, int] = {}
    septa_cache: Dict[str, int] = {
        s.septatype_string: s.septatype_id for s in session.scalars(select(Septatype)).all()
    }

    for r in df.to_dict(orient="records"):
        gid = r["Genome ID"]
        kwargs: Dict[str, Any] = dict(
            bv_brc_genome_id=gid,
            isolate_name=r["Isolate Name"],
            species_id=_get_or_create_species(session, sp_cache, r["species"]),
            habitat_id=_get_or_create_habitat(session, hab_cache, r["habitat"]),
            mlst_st=r["MLST"],
            collection_year_start=r["year_start"],
            collection_year_end=r["year_end"],
            collection_year_raw=r["year_raw"],
            country_name=r["Isolation Country"],
            country_iso3=r["country_iso3"],
            lab_phenotype=r["lab_phenotype"],
            septatype_id=_get_or_create_septatype(session, septa_cache, r),
            source=IsolateSource.FOUNDING,
        )
        if bvbrc and gid in bvbrc:
            rec = bvbrc[gid]
            if rec.get("genome_name"):
                kwargs["official_genome_name"] = rec["genome_name"]
                stats["name_ok"] += 1
            else:
                stats["name_missing"] += 1
            pred, raw, method = summarise_computational(rec.get("amr", []))
            kwargs.update(
                bv_brc_computational_prediction=pred,
                bv_brc_prediction_raw=raw,
                bv_brc_prediction_method=method,
                bv_brc_fetched_at=datetime.fromisoformat(rec["fetched_at"]) if rec.get("fetched_at") else None,
            )
            stats[f"bvbrc_{pred.value}"] += 1
        elif bvbrc is not None:
            stats["api_not_fetched"] += 1

        session.merge(Isolate(**kwargs))  # unset attributes are left untouched on re-runs
        session.flush()

        # Genotype: amino acids only; codon columns are never touched here.
        session.merge(IsolateGenotype(
            bv_brc_genome_id=gid,
            gyra83_aa=r["GyrA-83"], gyra87_aa=r["GyrA-87"], parc80_aa=r["ParC-80"],
            parc84_aa=r["ParC-84"], recn346_aa=r["RecN-346"],
        ))

        # PMQR calls: replace this isolate's qnrB2/qnrS1 rows
        session.query(PMQRCall).filter(
            PMQRCall.bv_brc_genome_id == gid, PMQRCall.gene.in_(["qnrB2", "qnrS1"])
        ).delete(synchronize_session=False)
        for gene in ("qnrB2", "qnrS1"):
            if r[gene] is not None:
                session.add(PMQRCall(bv_brc_genome_id=gid, gene=gene, present=bool(r[gene])))
        stats["isolates"] += 1

    session.flush()
    # refresh septatype isolate counts
    counts = (
        select(func.count(Isolate.bv_brc_genome_id))
        .where(Isolate.septatype_id == Septatype.septatype_id)
        .scalar_subquery()
    )
    session.execute(update(Septatype).values(isolate_count=counts))
    session.commit()
    return stats


# =============================================================================
# 4. CLI
# =============================================================================
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True, type=Path, help="Path to Table S1.xlsx")
    ap.add_argument("--database-url", default=None, help="Overrides $DATABASE_URL")
    ap.add_argument("--api-base", default="https://www.bv-brc.org/api")
    ap.add_argument("--sleep", type=float, default=0.5, help="Seconds between API calls (default 0.5)")
    ap.add_argument("--batch-size", type=int, default=50, help="Genome IDs per API request; use 1 for strictly one-by-one")
    ap.add_argument("--cache", type=Path, default=Path("bvbrc_cache.json"))
    ap.add_argument("--skip-api", action="store_true", help="Load Table S1 only, no BV-BRC calls")
    ap.add_argument("--limit", type=int, default=None, help="Only process the first N isolates (testing)")
    ap.add_argument("--dry-run", action="store_true", help="Parse and report only; touch nothing")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler("ingest.log")])

    df, warnings = parse_table_s1(args.xlsx)
    if args.limit:
        df = df.head(args.limit)
    log.info("Parsed %d isolates", len(df))
    log.info("By species/habitat:\n%s", df.groupby(["species", "habitat"]).size().to_string())
    log.info("Distinct septatypes: %d | null MLST: %d | phenotype: %s",
             df["septatype_string"].nunique(), df["MLST"].isna().sum(),
             df["lab_phenotype"].map(lambda p: p.value).value_counts().to_dict())
    for w in warnings:
        log.warning(w)
    Path("ingest_warnings.txt").write_text("\n".join(warnings))
    if args.dry_run:
        log.info("Dry run complete (%d warnings, written to ingest_warnings.txt)", len(warnings))
        return 0

    bvbrc = None
    if not args.skip_api:
        client = BVBRCClient(base=args.api_base, sleep=args.sleep, batch_size=args.batch_size)
        bvbrc = client.fetch_all(df["Genome ID"].tolist(), args.cache)

    init_db(args.database_url)
    get_engine(args.database_url)
    with SessionLocal() as session:
        seed_lookups(session)
        stats = load_isolates(session, df, bvbrc)

    log.info("Load complete: %s", dict(stats))
    if bvbrc is not None and stats.get("bvbrc_No_call", 0) == stats.get("isolates", 0):
        log.warning("Every isolate came back No_call. Check the antibiotic name and evidence values "
                    "in %s against a live genome_amr record.", args.cache)
    return 0


if __name__ == "__main__":
    sys.exit(main())
