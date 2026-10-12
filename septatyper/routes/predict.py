"""POST /api/v1/predict: upload a genome assembly, get a SeptaTyper prediction."""
import os
from typing import Dict, List, Literal, Optional

from fastapi import APIRouter, Depends, File, Form, HTTPException, Request, UploadFile
from pydantic import BaseModel
from sqlalchemy import select
from sqlalchemy.orm import Session

from ..database import get_session
from ..models import Septatype
from ..engine_runtime import EngineBusy, get_prediction_engine, run_engine  # noqa: F401  (re-exported)
from ..prediction_engine import EngineError, InputError
from ..uploads import read_fasta_upload

router = APIRouter(tags=["prediction"])

MAX_UPLOAD_BYTES = int(os.getenv("SEPTATYPER_MAX_UPLOAD_BYTES", str(60 * 1024 * 1024)))

SpeciesName = Literal[
    "Escherichia coli", "Shigella sonnei", "Salmonella enterica subspecies I", "Klebsiella pneumoniae"
]


# ------------------------------------------------------------------ response models
class CodonOut(BaseModel):
    position: str
    gene: str
    ref_nt_start: int
    codon: Optional[str] = None
    aa: Optional[str] = None
    status: str


class RecN346Out(BaseModel):
    codon: Optional[str] = None
    aa: Optional[str] = None
    leu_codon: Optional[str] = None
    synonymous_flag: Optional[bool] = None


class PMQROut(BaseModel):
    present: int
    identity: Optional[float] = None
    coverage: Optional[float] = None
    contig: Optional[str] = None
    closest_allele: Optional[str] = None


class SpeciesOut(BaseModel):
    name: Optional[str] = None
    source: str
    genes_agreeing: Optional[str] = None
    inferred: Optional[str] = None


class Tier1Out(BaseModel):
    baseline: str
    votes: int
    decided_by: str
    positions: List[dict]


class Tier2Out(BaseModel):
    override_applied: bool
    rules: List[str]


class DatabaseMatch(BaseModel):
    in_database: bool
    septatype_id: Optional[int] = None
    isolate_count: Optional[int] = None
    is_novel: bool


class PredictionResponse(BaseModel):
    pipeline_version: str
    input: Dict[str, int]
    species: SpeciesOut
    septatype: Optional[str] = None
    septatype_aa_level: Optional[str] = None
    prediction: Literal["Cres", "Csus", "Indeterminate"]
    predicted_phenotype: Optional[Literal["Resistant", "Susceptible"]] = None
    mean_hydrophobicity: Optional[float] = None
    hydrophobicity_scale: str
    codons: Dict[str, CodonOut]
    recn346: RecN346Out
    pmqr: Dict[str, PMQROut]
    genes: Dict[str, dict]
    tier1: Optional[Tier1Out] = None
    tier2: Optional[Tier2Out] = None
    database: Optional[DatabaseMatch] = None
    warnings: List[str]


# --------------------------------------------------------------------------- route
@router.post("/predict", response_model=PredictionResponse, summary="Predict ciprofloxacin phenotype from a genome")
def predict(
    request: Request,
    file: Optional[UploadFile] = File(None, description="Genome assembly, nucleotide FASTA (.fa/.fasta, optionally .gz)"),
    sequence: Optional[str] = Form(None, description="Alternatively, paste FASTA text"),
    species: Optional[SpeciesName] = Form(None, description="Recommended. Needed to apply the S. sonnei rule reliably"),
    db: Session = Depends(get_session),
):
    declared = request.headers.get("content-length")
    if declared and declared.isdigit() and int(declared) > MAX_UPLOAD_BYTES + 1024 * 1024:
        raise HTTPException(413, "Request body too large")
    if (file is None) == (not sequence):
        raise HTTPException(422, "Provide exactly one of 'file' or 'sequence'")

    text = read_fasta_upload(file, MAX_UPLOAD_BYTES) if file is not None else (sequence or "")
    if len(text) > MAX_UPLOAD_BYTES:
        raise HTTPException(413, "Sequence text too large")
    if not text.lstrip().startswith(">"):
        text = ">pasted\n" + text

    try:
        get_prediction_engine()
    except EngineError as exc:
        raise HTTPException(503, f"Prediction engine unavailable: {exc}") from exc

    try:
        result = run_engine(text, species)
    except EngineBusy as exc:
        raise HTTPException(503, "The server is busy. Try again in a minute.") from exc
    except InputError as exc:
        raise HTTPException(400, str(exc)) from exc
    except EngineError as exc:
        raise HTTPException(500, f"Pipeline error: {exc}") from exc

    st = result["septatype"]
    if st:
        row = db.scalar(select(Septatype).where(Septatype.septatype_string == st))
        result["database"] = {
            "in_database": row is not None,
            "septatype_id": row.septatype_id if row else None,
            "isolate_count": row.isolate_count if row else None,
            "is_novel": row is None,
        }
    return result
