"""Pydantic response models for the public API."""
from __future__ import annotations

import re
from typing import List, Optional

from pydantic import BaseModel

BVBRC_GENOME_URL = "https://www.bv-brc.org/view/Genome/"
BVBRC_ID_RE = re.compile(r"^\d+\.\d+$")


def bvbrc_url(genome_id: str) -> Optional[str]:
    return BVBRC_GENOME_URL + genome_id if BVBRC_ID_RE.match(genome_id) else None


class IsolateRow(BaseModel):
    bv_brc_genome_id: str
    bv_brc_url: Optional[str] = None  # None for community isolates without a BV-BRC ID
    official_genome_name: Optional[str] = None
    isolate_name: Optional[str] = None
    species: str
    habitat: Optional[str] = None
    mlst_st: Optional[str] = None
    septatype: Optional[str] = None
    collection_year_raw: Optional[str] = None
    collection_year_start: Optional[int] = None
    collection_year_end: Optional[int] = None
    country_name: Optional[str] = None
    country_iso3: Optional[str] = None
    lab_phenotype: str
    bv_brc_computational_prediction: Optional[str] = None
    septatyper_prediction: Optional[str] = None
    septatyper_correct: Optional[bool] = None
    bvbrc_correct: Optional[bool] = None
    septatyper_wins: Optional[bool] = None


class IsolatePage(BaseModel):
    total: int
    page: int
    page_size: int
    items: List[IsolateRow]


class PMQRItem(BaseModel):
    gene: str
    present: bool
    identity: Optional[float] = None
    coverage: Optional[float] = None


class GenotypeOut(BaseModel):
    gyra83_codon: Optional[str] = None
    gyra83_aa: Optional[str] = None
    gyra87_codon: Optional[str] = None
    gyra87_aa: Optional[str] = None
    parc80_codon: Optional[str] = None
    parc80_aa: Optional[str] = None
    parc84_codon: Optional[str] = None
    parc84_aa: Optional[str] = None
    recn346_codon: Optional[str] = None
    recn346_aa: Optional[str] = None
    recn346_synonymous_flag: Optional[bool] = None


class IsolateDetail(IsolateRow):
    bv_brc_prediction_raw: Optional[str] = None
    bv_brc_prediction_method: Optional[str] = None
    mean_hydrophobicity: Optional[float] = None
    genotype: Optional[GenotypeOut] = None
    pmqr: List[PMQRItem] = []


class FacetItem(BaseModel):
    value: str
    count: int


class Facets(BaseModel):
    species: List[FacetItem]
    habitat: List[FacetItem]
    st: List[FacetItem]
    septatype: List[FacetItem]
    country: List[FacetItem]
    year_min: Optional[int] = None
    year_max: Optional[int] = None


class Summary(BaseModel):
    total: int
    lab_resistant: int
    lab_susceptible: int
    septatyper_predicted: int
    septatyper_accuracy: Optional[float] = None
    bvbrc_predicted: int
    bvbrc_accuracy: Optional[float] = None
    bvbrc_no_call: int
    septatyper_wins: int


class CountryStat(BaseModel):
    iso3: str
    country: Optional[str] = None
    n_isolates: int
    n_resistant: int
    resistant_fraction: float


class HabitatSeptatypeStat(BaseModel):
    habitat: str
    septatype: str
    n: int
