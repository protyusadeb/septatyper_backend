from typing import Annotated, Literal

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from ..database import get_session
from ..models import (
    Habitat,
    Isolate,
    IsolateGenotype,
    PMQRCall,
    Septatype,
    Species,
)
from ..schemas import (
    bvbrc_url,
    Facets,
    FacetItem,
    GenotypeOut,
    IsolateDetail,
    IsolatePage,
    IsolateRow,
    PMQRItem,
)
from .common import IsolateFilters, apply_filters, base_select

router = APIRouter(tags=["isolates"])

SortKey = Literal[
    "bv_brc_genome_id", "species", "habitat", "mlst_st", "septatype", "collection_year",
    "country", "lab_phenotype", "bv_brc_prediction", "septatyper_prediction",
]
SORT_COLUMNS = {
    "bv_brc_genome_id": Isolate.bv_brc_genome_id,
    "species": Species.scientific_name,
    "habitat": Habitat.name,
    "mlst_st": Isolate.mlst_st,
    "septatype": Septatype.septatype_string,
    "collection_year": func.coalesce(Isolate.collection_year_end, Isolate.collection_year_start),
    "country": Isolate.country_name,
    "lab_phenotype": Isolate.lab_phenotype,
    "bv_brc_prediction": Isolate.bv_brc_computational_prediction,
    "septatyper_prediction": Isolate.septatyper_prediction,
}

ROW_COLUMNS = (
    Isolate.bv_brc_genome_id,
    Isolate.official_genome_name,
    Isolate.isolate_name,
    Species.scientific_name.label("species"),
    Habitat.name.label("habitat"),
    Isolate.mlst_st,
    Septatype.septatype_string.label("septatype"),
    Isolate.collection_year_raw,
    Isolate.collection_year_start,
    Isolate.collection_year_end,
    Isolate.country_name,
    Isolate.country_iso3,
    Isolate.lab_phenotype,
    Isolate.bv_brc_computational_prediction,
    Isolate.septatyper_prediction,
    Isolate.septatyper_correct,
    Isolate.bvbrc_correct,
    Isolate.septatyper_wins,
)


def _val(x):
    """Enum member -> plain string; everything else unchanged."""
    return x.value if hasattr(x, "value") else x


def row_to_schema(m) -> IsolateRow:
    d = {k: _val(v) for k, v in m.items()}
    d["bv_brc_url"] = bvbrc_url(d["bv_brc_genome_id"])
    return IsolateRow(**d)


@router.get("/isolates", response_model=IsolatePage, summary="Search and filter isolates")
def list_isolates(
    filters: Annotated[IsolateFilters, Depends()],
    page: Annotated[int, Query(ge=1)] = 1,
    page_size: Annotated[int, Query(ge=1, le=200)] = 25,
    sort_by: SortKey = "bv_brc_genome_id",
    sort_dir: Literal["asc", "desc"] = "asc",
    db: Session = Depends(get_session),
):
    count_stmt = apply_filters(base_select(func.count(Isolate.bv_brc_genome_id)), filters)
    total = db.scalar(count_stmt) or 0

    col = SORT_COLUMNS[sort_by]
    order = col.desc().nulls_last() if sort_dir == "desc" else col.asc().nulls_last()
    stmt = (
        apply_filters(base_select(*ROW_COLUMNS), filters)
        .order_by(order, Isolate.bv_brc_genome_id)  # PK tiebreak keeps paging stable
        .limit(page_size)
        .offset((page - 1) * page_size)
    )
    rows = db.execute(stmt).mappings().all()
    return IsolatePage(total=total, page=page, page_size=page_size, items=[row_to_schema(r) for r in rows])


@router.get("/isolates/{genome_id}", response_model=IsolateDetail, summary="Single isolate with genotype")
def get_isolate(genome_id: str, db: Session = Depends(get_session)):
    stmt = base_select(
        *ROW_COLUMNS,
        Isolate.bv_brc_prediction_raw,
        Isolate.bv_brc_prediction_method,
        Septatype.mean_hydrophobicity,
    ).where(Isolate.bv_brc_genome_id == genome_id)
    m = db.execute(stmt).mappings().first()
    if m is None:
        raise HTTPException(status_code=404, detail=f"Genome ID {genome_id} not found")
    d = {k: _val(v) for k, v in m.items()}
    d["bv_brc_url"] = bvbrc_url(genome_id)
    geno = db.get(IsolateGenotype, genome_id)
    d["genotype"] = (
        GenotypeOut(**{c: getattr(geno, c) for c in GenotypeOut.model_fields}) if geno else None
    )
    calls = db.scalars(select(PMQRCall).where(PMQRCall.bv_brc_genome_id == genome_id)).all()
    d["pmqr"] = [PMQRItem(gene=c.gene, present=c.present, identity=c.identity, coverage=c.coverage) for c in calls]
    return IsolateDetail(**d)


@router.get("/filters", response_model=Facets, summary="Filter options with counts")
def get_filters(db: Session = Depends(get_session)):
    def facet(col) -> list[FacetItem]:
        rows = db.execute(
            base_select(col, func.count()).where(col.is_not(None)).group_by(col).order_by(func.count().desc(), col)
        ).all()
        return [FacetItem(value=str(v), count=n) for v, n in rows]

    years = db.execute(
        select(
            func.min(func.coalesce(Isolate.collection_year_start, Isolate.collection_year_end)),
            func.max(func.coalesce(Isolate.collection_year_end, Isolate.collection_year_start)),
        )
    ).one()
    return Facets(
        species=facet(Species.scientific_name),
        habitat=facet(Habitat.name),
        st=facet(Isolate.mlst_st),
        septatype=facet(Septatype.septatype_string),
        country=facet(Isolate.country_iso3),
        year_min=years[0],
        year_max=years[1],
    )
