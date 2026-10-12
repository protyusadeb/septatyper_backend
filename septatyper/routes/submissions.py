"""Module C: community submissions (public) and the curation queue (curators/admins)."""
import gzip
import hashlib
import hmac
import os
import re
import secrets
import uuid
from datetime import date
from pathlib import Path
from typing import Annotated, Any, Dict, List, Literal, Optional

from fastapi import APIRouter, BackgroundTasks, Depends, File, Form, HTTPException, Query, Request, UploadFile
from pydantic import BaseModel, Field
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from ..auth import require_curator
from ..countries import country_to_iso3
from ..database import get_session
from ..engine_runtime import get_prediction_engine
from ..models import Isolate, Phenotype, Submission, SubmissionStatus, User
from ..prediction_engine import EngineError
from ..ratelimit import RateLimiter, enforce
from ..schemas import BVBRC_ID_RE
from ..submission_service import (
    ApprovalError, approve_submission, flag_objects, process_submission, reject_submission,
)
from ..uploads import looks_like_nucleotide_fasta, read_fasta_upload

router = APIRouter(tags=["submissions"])

MAX_UPLOAD_BYTES = int(os.getenv("SEPTATYPER_MAX_UPLOAD_BYTES", str(60 * 1024 * 1024)))
SUBMISSION_DIR = Path(os.getenv("SEPTATYPER_SUBMISSION_DIR", "submissions_data"))
_TRUST_PROXY = os.getenv("SEPTATYPER_TRUST_PROXY") == "1"
submit_limiter = RateLimiter(int(os.getenv("SEPTATYPER_SUBMIT_LIMIT", "10")), 3600)

SpeciesName = Literal[
    "Escherichia coli", "Shigella sonnei", "Salmonella enterica subspecies I", "Klebsiella pneumoniae"
]
HabitatName = Literal["Fecal", "Urine", "Respiratory", "Other"]
PhenotypeName = Literal["Resistant", "Susceptible"]
EMAIL_RE = re.compile(r"^[^@\s]{1,64}@[^@\s]{1,255}\.[^@\s]{2,}$")
_CTRL = re.compile(r"[\x00-\x1f\x7f]")


# ------------------------------------------------------------------------- models
class FlagOut(BaseModel):
    code: str
    label: str
    severity: str


class SubmissionAccepted(BaseModel):
    submission_id: int
    status_token: str = Field(description="Shown once. Needed to check progress later.")
    analysis_state: str


class SubmissionStatusOut(BaseModel):
    submission_id: int
    status: str
    analysis_state: Optional[str] = None
    analysis_error: Optional[str] = None
    bv_brc_id_supplied: Optional[str] = None
    bv_brc_lookup_status: Optional[str] = None
    bv_brc_name_fetched: Optional[str] = None
    septatype: Optional[str] = None
    septatyper_prediction: Optional[str] = None
    experimental_phenotype: Optional[str] = None
    flags: List[FlagOut] = []


class QueueItem(BaseModel):
    submission_id: int
    created_at: Optional[str] = None
    status: str
    analysis_state: Optional[str] = None
    isolate_name: Optional[str] = None
    declared_species: Optional[str] = None
    habitat: Optional[str] = None
    country: Optional[str] = None
    collection_year: Optional[int] = None
    experimental_phenotype: Optional[str] = None
    septatyper_prediction: Optional[str] = None
    predicted_septatype_string: Optional[str] = None
    is_novel_septatype: bool = False
    prediction_matches_phenotype: Optional[bool] = None
    bv_brc_id_supplied: Optional[str] = None
    bv_brc_name_fetched: Optional[str] = None
    flags: List[FlagOut] = []


class QueueCounts(BaseModel):
    needs_review: int
    pending: int
    failed: int


class QueuePage(BaseModel):
    total: int
    page: int
    page_size: int
    counts: QueueCounts
    items: List[QueueItem]


class SubmissionDetail(QueueItem):
    contact_email: Optional[str] = None
    submitter_name: Optional[str] = None
    source: Optional[str] = None
    experimental_mic: Optional[float] = None
    mic_unit: Optional[str] = None
    mic_method: Optional[str] = None
    bv_brc_lookup_status: Optional[str] = None
    analysis_error: Optional[str] = None
    prediction: Optional[Dict[str, Any]] = None
    curator_notes: Optional[str] = None
    decided_at: Optional[str] = None


class ApproveRequest(BaseModel):
    notes: Optional[str] = Field(None, max_length=2000)
    confirm_flags: bool = False


class RejectRequest(BaseModel):
    notes: str = Field(min_length=3, max_length=2000)


class ApproveResponse(BaseModel):
    isolate_id: str
    septatype_id: int
    septatype_string: str
    created_septatype: bool


# ------------------------------------------------------------------------ helpers
def _text(value: Optional[str], name: str, max_len: int) -> Optional[str]:
    if value is None:
        return None
    v = _CTRL.sub(" ", value).strip()
    if not v:
        return None
    if len(v) > max_len:
        raise HTTPException(422, f"{name} is too long (max {max_len} characters)")
    return v


def _qitem(s: Submission) -> dict:
    return dict(
        submission_id=s.submission_id,
        created_at=s.created_at.isoformat() if s.created_at else None,
        status=s.status.value,
        analysis_state=s.analysis_state,
        isolate_name=s.isolate_name,
        declared_species=s.declared_species,
        habitat=s.habitat,
        country=s.country,
        collection_year=s.collection_year,
        experimental_phenotype=s.experimental_phenotype.value if s.experimental_phenotype else None,
        septatyper_prediction=s.septatyper_prediction,
        predicted_septatype_string=s.predicted_septatype_string,
        is_novel_septatype=bool(s.is_novel_septatype),
        prediction_matches_phenotype=s.prediction_matches_phenotype,
        bv_brc_id_supplied=s.bv_brc_id_supplied,
        bv_brc_name_fetched=s.bv_brc_name_fetched,
        flags=flag_objects(s.flags),
    )


def _sequence_hash(text: str) -> str:
    """Hash of the bases only, so re-wrapped or renamed copies of one genome still collide."""
    seq = "".join(l.strip() for l in text.splitlines() if l and not l.startswith(">")).upper()
    return hashlib.sha256(seq.encode()).hexdigest()


# ------------------------------------------------------------------------- public
@router.post("/submissions", response_model=SubmissionAccepted, status_code=202, summary="Submit a genome with metadata")
def create_submission(
    request: Request,
    background_tasks: BackgroundTasks,
    file: Annotated[UploadFile, File(description="Genome assembly, nucleotide FASTA (optionally .gz)")],
    contact_email: Annotated[str, Form()],
    species: Annotated[SpeciesName, Form()],
    habitat: Annotated[HabitatName, Form()],
    country: Annotated[str, Form()],
    collection_year: Annotated[int, Form()],
    lab_phenotype: Annotated[PhenotypeName, Form(description="Ciprofloxacin phenotype from the lab")],
    submitter_name: Annotated[Optional[str], Form()] = None,
    mic_value: Annotated[Optional[str], Form(description="Ciprofloxacin MIC")] = None,
    mic_unit: Annotated[Optional[str], Form()] = None,
    mic_method: Annotated[Optional[str], Form()] = None,
    bv_brc_id: Annotated[Optional[str], Form(description="Optional, e.g. 562.7710")] = None,
    isolate_name: Annotated[Optional[str], Form()] = None,
    source: Annotated[Optional[str], Form()] = None,
    db: Session = Depends(get_session),
):
    enforce(submit_limiter, request, _TRUST_PROXY)

    email = (contact_email or "").strip().lower()
    if not EMAIL_RE.match(email):
        raise HTTPException(422, "Enter a valid contact email address")
    if not 1900 <= collection_year <= date.today().year:
        raise HTTPException(422, f"Collection year must be between 1900 and {date.today().year}")
    country_clean = _text(country, "Country", 80)
    iso3 = country_to_iso3(country_clean) if country_clean else None
    if iso3 is None:
        raise HTTPException(422, "Country not recognised. Use the full country name or its 3-letter ISO code")

    mic: Optional[float] = None
    unit = _text(mic_unit, "MIC unit", 10)
    mv = _text(mic_value, "MIC", 20)
    if mv is not None:
        try:
            mic = float(mv)
        except ValueError:
            raise HTTPException(422, "MIC must be a number")
        if not 0 < mic <= 1024:
            raise HTTPException(422, "MIC must be greater than 0 and at most 1024")
        unit = unit or "mg/L"
        if unit not in ("mg/L", "ug/mL"):
            raise HTTPException(422, "MIC unit must be mg/L or ug/mL")
    else:
        unit = None

    bvid = _text(bv_brc_id, "BV-BRC ID", 30)
    if bvid is not None and not BVBRC_ID_RE.match(bvid):
        raise HTTPException(422, "BV-BRC Genome IDs look like 562.7710")

    text = read_fasta_upload(file, MAX_UPLOAD_BYTES)
    if not looks_like_nucleotide_fasta(text):
        raise HTTPException(400, "Upload a nucleotide FASTA file (it must start with a '>' header line)")
    try:
        get_prediction_engine()
    except EngineError as exc:
        raise HTTPException(503, f"Analysis is temporarily unavailable: {exc}") from exc

    if bvid:
        if db.get(Isolate, bvid) is not None:
            raise HTTPException(409, f"Genome {bvid} is already in the database")
        if db.scalar(select(func.count()).select_from(Submission).where(
            Submission.bv_brc_id_supplied == bvid,
            Submission.status.in_([SubmissionStatus.PENDING, SubmissionStatus.NEEDS_REVIEW, SubmissionStatus.APPROVED]),
        )):
            raise HTTPException(409, f"Genome {bvid} has already been submitted")

    SUBMISSION_DIR.mkdir(parents=True, exist_ok=True, mode=0o700)
    path = SUBMISSION_DIR / f"{uuid.uuid4().hex}.fa.gz"  # never derived from the uploaded filename
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(text)
    os.chmod(path, 0o600)

    token = secrets.token_urlsafe(24)
    sub = Submission(
        status=SubmissionStatus.PENDING, analysis_state="queued",
        contact_email=email, submitter_name=_text(submitter_name, "Name", 120),
        isolate_name=_text(isolate_name, "Isolate name", 200), source=_text(source, "Source", 80),
        declared_species=species, habitat=habitat, country=country_clean, country_iso3=iso3,
        collection_year=collection_year, experimental_phenotype=Phenotype(lab_phenotype),
        experimental_mic=mic, mic_unit=unit, mic_method=_text(mic_method, "MIC method", 60),
        bv_brc_id_supplied=bvid, file_path=str(path), file_sha256=_sequence_hash(text),
        status_token_hash=hashlib.sha256(token.encode()).hexdigest(),
    )
    db.add(sub)
    db.commit()
    background_tasks.add_task(process_submission, db.get_bind(), sub.submission_id)
    return SubmissionAccepted(submission_id=sub.submission_id, status_token=token, analysis_state="queued")


@router.get("/submissions/{submission_id}/status", response_model=SubmissionStatusOut,
            summary="Submitter progress check (needs the token from the upload response)")
def submission_status(submission_id: int, token: Annotated[str, Query(min_length=10, max_length=100)],
                      db: Session = Depends(get_session)):
    sub = db.get(Submission, submission_id)
    given = hashlib.sha256(token.encode()).hexdigest()
    if sub is None or not sub.status_token_hash or not hmac.compare_digest(given, sub.status_token_hash):
        raise HTTPException(404, "Submission not found")
    return SubmissionStatusOut(
        submission_id=sub.submission_id, status=sub.status.value, analysis_state=sub.analysis_state,
        analysis_error=sub.analysis_error, bv_brc_id_supplied=sub.bv_brc_id_supplied,
        bv_brc_lookup_status=sub.bv_brc_lookup_status, bv_brc_name_fetched=sub.bv_brc_name_fetched,
        septatype=sub.predicted_septatype_string, septatyper_prediction=sub.septatyper_prediction,
        experimental_phenotype=sub.experimental_phenotype.value if sub.experimental_phenotype else None,
        flags=flag_objects(sub.flags),
    )


# ------------------------------------------------------------------------ curators
@router.get("/submissions/pending", response_model=QueuePage, summary="Curation queue (curators and admins)")
def pending_queue(
    status: Literal["all", "needs_review", "pending"] = "all",
    novel: Optional[bool] = None,
    mismatch: Optional[bool] = Query(None, description="true = only phenotype mismatches"),
    page: Annotated[int, Query(ge=1)] = 1,
    page_size: Annotated[int, Query(ge=1, le=100)] = 25,
    _: User = Depends(require_curator),
    db: Session = Depends(get_session),
):
    open_states = [SubmissionStatus.PENDING, SubmissionStatus.NEEDS_REVIEW]
    base = select(Submission).where(
        Submission.status.in_(open_states), Submission.analysis_state.in_(["done", "failed"])
    )
    stmt = base
    if status != "all":
        stmt = stmt.where(Submission.status == SubmissionStatus(status))
    if novel is not None:
        stmt = stmt.where(Submission.is_novel_septatype.is_(novel))
    if mismatch:
        stmt = stmt.where(Submission.prediction_matches_phenotype.is_(False))
    total = db.scalar(select(func.count()).select_from(stmt.subquery())) or 0
    rows = db.scalars(
        stmt.order_by(Submission.created_at.asc(), Submission.submission_id.asc())
        .limit(page_size).offset((page - 1) * page_size)
    ).all()

    def count(cond) -> int:
        return db.scalar(select(func.count()).select_from(base.where(cond).subquery())) or 0

    counts = QueueCounts(
        needs_review=count(Submission.status == SubmissionStatus.NEEDS_REVIEW),
        pending=count(Submission.status == SubmissionStatus.PENDING),
        failed=count(Submission.analysis_state == "failed"),
    )
    return QueuePage(total=total, page=page, page_size=page_size, counts=counts, items=[QueueItem(**_qitem(s)) for s in rows])


@router.get("/submissions/{submission_id}", response_model=SubmissionDetail, summary="Full submission detail (curators)")
def submission_detail(submission_id: int, _: User = Depends(require_curator), db: Session = Depends(get_session)):
    sub = db.get(Submission, submission_id)
    if sub is None:
        raise HTTPException(404, "Submission not found")
    d = _qitem(sub)
    d.update(
        contact_email=sub.contact_email, submitter_name=sub.submitter_name, source=sub.source,
        experimental_mic=float(sub.experimental_mic) if sub.experimental_mic is not None else None,
        mic_unit=sub.mic_unit, mic_method=sub.mic_method, bv_brc_lookup_status=sub.bv_brc_lookup_status,
        analysis_error=sub.analysis_error, prediction=sub.prediction_result, curator_notes=sub.curator_notes,
        decided_at=sub.decided_at.isoformat() if sub.decided_at else None,
    )
    return SubmissionDetail(**d)


@router.post("/submissions/{submission_id}/approve", response_model=ApproveResponse,
             summary="Approve and copy into the live database (curators)")
def approve(submission_id: int, body: ApproveRequest, curator: User = Depends(require_curator),
            db: Session = Depends(get_session)):
    try:
        out = approve_submission(db, submission_id, curator, (body.notes or "").strip() or None, body.confirm_flags)
        db.commit()
    except ApprovalError as exc:
        db.rollback()
        raise HTTPException(exc.status_code, exc.message) from exc
    return ApproveResponse(**out)


@router.post("/submissions/{submission_id}/reject", status_code=204, summary="Reject a submission (curators)")
def reject(submission_id: int, body: RejectRequest, curator: User = Depends(require_curator),
           db: Session = Depends(get_session)):
    try:
        reject_submission(db, submission_id, curator, body.notes.strip())
        db.commit()
    except ApprovalError as exc:
        db.rollback()
        raise HTTPException(exc.status_code, exc.message) from exc
