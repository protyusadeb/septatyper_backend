"""Shared filter dependency and query builder used by the table and the charts."""
from typing import Annotated, List, Literal, Optional

from fastapi import Query
from sqlalchemy import Select, func, or_, select

from ..models import (
    BvbrcPrediction,
    Habitat,
    Isolate,
    Phenotype,
    SeptaPrediction,
    Septatype,
    Species,
)


class IsolateFilters:
    """Query parameters accepted by /isolates and every /stats endpoint.

    List parameters are repeated: ?species=Escherichia%20coli&species=Shigella%20sonnei
    """

    def __init__(
        self,
        q: Annotated[Optional[str], Query(max_length=100, description="Genome ID, genome name or isolate name")] = None,
        species: Annotated[Optional[List[str]], Query()] = None,
        habitat: Annotated[Optional[List[str]], Query()] = None,
        st: Annotated[Optional[List[str]], Query(description="Sequence type, e.g. ST131")] = None,
        septatype: Annotated[Optional[List[str]], Query(description="e.g. S-D-S-E-L-0-0")] = None,
        country: Annotated[Optional[List[str]], Query(description="ISO3 codes")] = None,
        year_min: Annotated[Optional[int], Query(ge=1900, le=2100)] = None,
        year_max: Annotated[Optional[int], Query(ge=1900, le=2100)] = None,
        lab_phenotype: Optional[Literal["Resistant", "Susceptible"]] = None,
        bvbrc_prediction: Optional[Literal["Resistant", "Susceptible", "No_call"]] = None,
        septatyper_prediction: Optional[Literal["Cres", "Csus"]] = None,
        septatyper_wins: Annotated[Optional[bool], Query(description="SeptaTyper right, BV-BRC wrong or no call")] = None,
        septatyper_correct: Optional[bool] = None,
        bvbrc_correct: Optional[bool] = None,
    ):
        self.q = q.strip() if q else None
        self.species, self.habitat, self.st = species, habitat, st
        self.septatype, self.country = septatype, country
        self.year_min, self.year_max = year_min, year_max
        self.lab_phenotype = lab_phenotype
        self.bvbrc_prediction = bvbrc_prediction
        self.septatyper_prediction = septatyper_prediction
        self.septatyper_wins = septatyper_wins
        self.septatyper_correct = septatyper_correct
        self.bvbrc_correct = bvbrc_correct


def base_select(*columns) -> Select:
    """SELECT <columns> FROM isolates JOIN species LEFT JOIN habitats LEFT JOIN septatypes."""
    return (
        select(*columns)
        .select_from(Isolate)
        .join(Species, Isolate.species_id == Species.species_id)
        .outerjoin(Habitat, Isolate.habitat_id == Habitat.habitat_id)
        .outerjoin(Septatype, Isolate.septatype_id == Septatype.septatype_id)
    )


def _escape_like(s: str) -> str:
    return s.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def apply_filters(stmt: Select, f: IsolateFilters) -> Select:
    if f.q:
        like = f"%{_escape_like(f.q)}%"
        stmt = stmt.where(
            or_(
                Isolate.bv_brc_genome_id.like(f"{_escape_like(f.q)}%", escape="\\"),
                Isolate.official_genome_name.ilike(like, escape="\\"),
                Isolate.isolate_name.ilike(like, escape="\\"),
            )
        )
    if f.species:
        stmt = stmt.where(Species.scientific_name.in_(f.species))
    if f.habitat:
        stmt = stmt.where(Habitat.name.in_(f.habitat))
    if f.st:
        stmt = stmt.where(Isolate.mlst_st.in_(f.st))
    if f.septatype:
        stmt = stmt.where(Septatype.septatype_string.in_(f.septatype))
    if f.country:
        stmt = stmt.where(Isolate.country_iso3.in_([c.upper() for c in f.country]))
    # Year filter uses range overlap; isolates with no usable year drop out when it is set.
    start = func.coalesce(Isolate.collection_year_start, Isolate.collection_year_end)
    end = func.coalesce(Isolate.collection_year_end, Isolate.collection_year_start)
    if f.year_min is not None:
        stmt = stmt.where(end >= f.year_min)
    if f.year_max is not None:
        stmt = stmt.where(start <= f.year_max)
    if f.lab_phenotype:
        stmt = stmt.where(Isolate.lab_phenotype == Phenotype(f.lab_phenotype))
    if f.bvbrc_prediction:
        stmt = stmt.where(Isolate.bv_brc_computational_prediction == BvbrcPrediction(f.bvbrc_prediction))
    if f.septatyper_prediction:
        stmt = stmt.where(Isolate.septatyper_prediction == SeptaPrediction(f.septatyper_prediction))
    for flag, col in (
        (f.septatyper_wins, Isolate.septatyper_wins),
        (f.septatyper_correct, Isolate.septatyper_correct),
        (f.bvbrc_correct, Isolate.bvbrc_correct),
    ):
        if flag is not None:
            stmt = stmt.where(col.is_(flag))
    return stmt
