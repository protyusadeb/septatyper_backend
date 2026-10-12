from collections import defaultdict
from typing import Annotated, List

from fastapi import APIRouter, Depends, Query
from sqlalchemy import case, func
from sqlalchemy.orm import Session

from ..database import get_session
from ..models import BvbrcPrediction, Habitat, Isolate, Phenotype, Septatype
from ..schemas import CountryStat, HabitatSeptatypeStat, Summary
from .common import IsolateFilters, apply_filters, base_select

router = APIRouter(prefix="/stats", tags=["stats"])


def _count_if(cond):
    return func.coalesce(func.sum(case((cond, 1), else_=0)), 0)


@router.get("/summary", response_model=Summary, summary="Counts and benchmark accuracy for the current filter")
def summary(filters: Annotated[IsolateFilters, Depends()], db: Session = Depends(get_session)):
    stmt = apply_filters(
        base_select(
            func.count(Isolate.bv_brc_genome_id),
            _count_if(Isolate.lab_phenotype == Phenotype.RESISTANT),
            _count_if(Isolate.lab_phenotype == Phenotype.SUSCEPTIBLE),
            _count_if(Isolate.septatyper_correct.is_not(None)),
            _count_if(Isolate.septatyper_correct.is_(True)),
            _count_if(Isolate.bvbrc_correct.is_not(None)),
            _count_if(Isolate.bvbrc_correct.is_(True)),
            _count_if(Isolate.bv_brc_computational_prediction == BvbrcPrediction.NO_CALL),
            _count_if(Isolate.septatyper_wins.is_(True)),
        ),
        filters,
    )
    total, res, sus, s_n, s_ok, b_n, b_ok, nocall, wins = db.execute(stmt).one()
    return Summary(
        total=total, lab_resistant=res, lab_susceptible=sus,
        septatyper_predicted=s_n, septatyper_accuracy=(s_ok / s_n) if s_n else None,
        bvbrc_predicted=b_n, bvbrc_accuracy=(b_ok / b_n) if b_n else None,
        bvbrc_no_call=nocall, septatyper_wins=wins,
    )


@router.get("/countries", response_model=List[CountryStat], summary="Per-country counts for the choropleth")
def countries(filters: Annotated[IsolateFilters, Depends()], db: Session = Depends(get_session)):
    stmt = (
        apply_filters(
            base_select(
                Isolate.country_iso3,
                func.min(Isolate.country_name),
                func.count(Isolate.bv_brc_genome_id),
                _count_if(Isolate.lab_phenotype == Phenotype.RESISTANT),
            ),
            filters,
        )
        .where(Isolate.country_iso3.is_not(None))
        .group_by(Isolate.country_iso3)
        .order_by(func.count(Isolate.bv_brc_genome_id).desc())
    )
    return [
        CountryStat(iso3=iso, country=name, n_isolates=n, n_resistant=r, resistant_fraction=r / n)
        for iso, name, n, r in db.execute(stmt).all()
    ]


@router.get("/septatype-habitat", response_model=List[HabitatSeptatypeStat],
            summary="Septatype counts per habitat (top N, rest grouped as 'Other')")
def septatype_by_habitat(
    filters: Annotated[IsolateFilters, Depends()],
    top: Annotated[int, Query(ge=1, le=20)] = 8,
    db: Session = Depends(get_session),
):
    stmt = (
        apply_filters(base_select(Habitat.name, Septatype.septatype_string, func.count()), filters)
        .where(Habitat.name.is_not(None), Septatype.septatype_string.is_not(None))
        .group_by(Habitat.name, Septatype.septatype_string)
    )
    rows = db.execute(stmt).all()
    overall: dict[str, int] = defaultdict(int)
    for _, s, n in rows:
        overall[s] += n
    keep = {s for s, _ in sorted(overall.items(), key=lambda kv: (-kv[1], kv[0]))[:top]}
    merged: dict[tuple[str, str], int] = defaultdict(int)
    for h, s, n in rows:
        merged[(h, s if s in keep else "Other")] += n
    return [HabitatSeptatypeStat(habitat=h, septatype=s, n=n) for (h, s), n in sorted(merged.items())]
