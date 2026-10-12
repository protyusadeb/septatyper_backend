"""SQLAlchemy 2.0 ORM models for SeptaTyper.

Design notes
- bv_brc_genome_id (e.g. "562.7710") is the primary key of isolates.
- Enums are stored as VARCHAR + CHECK constraint (native_enum=False). This keeps
  the generated benchmarking columns immutable, which PostgreSQL requires.
- Septatype position order everywhere:
  GyrA-83, GyrA-87, ParC-80, ParC-84, RecN-346, qnrB2, qnrS1.
- Codons are stored as uppercase DNA triplets (e.g. "CTA", "CTG").
"""
from __future__ import annotations

import datetime as dt
import enum
from typing import List, Optional

from sqlalchemy import (
    JSON,
    Boolean,
    CheckConstraint,
    Computed,
    Date,
    DateTime,
    Enum as SAEnum,
    Float,
    ForeignKey,
    Index,
    Integer,
    MetaData,
    Numeric,
    SmallInteger,
    String,
    Text,
    UniqueConstraint,
    func,
)
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship

JSONType = JSON().with_variant(JSONB(), "postgresql")


class Base(DeclarativeBase):
    metadata = MetaData(
        naming_convention={
            "ix": "ix_%(column_0_label)s",
            "uq": "uq_%(table_name)s_%(column_0_name)s",
            "ck": "ck_%(table_name)s_%(constraint_name)s",
            "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
            "pk": "pk_%(table_name)s",
        }
    )


# --------------------------------------------------------------------------- enums
class Phenotype(str, enum.Enum):
    RESISTANT = "Resistant"
    SUSCEPTIBLE = "Susceptible"


class BvbrcPrediction(str, enum.Enum):
    RESISTANT = "Resistant"
    SUSCEPTIBLE = "Susceptible"
    NO_CALL = "No_call"


class SeptaPrediction(str, enum.Enum):
    CRES = "Cres"
    CSUS = "Csus"


class IsolateSource(str, enum.Enum):
    FOUNDING = "founding"
    COMMUNITY = "community"


class SubmissionStatus(str, enum.Enum):
    PENDING = "pending"
    NEEDS_REVIEW = "needs_review"
    APPROVED = "approved"
    REJECTED = "rejected"


class UserRole(str, enum.Enum):
    PUBLIC = "public"
    SUBMITTER = "submitter"
    CURATOR = "curator"
    ADMIN = "admin"


class JobInputType(str, enum.Enum):
    FASTA = "fasta"
    READS = "reads"
    PASTED = "pasted"


class JobStatus(str, enum.Enum):
    QUEUED = "queued"
    RUNNING = "running"
    DONE = "done"
    FAILED = "failed"


class EffectType(str, enum.Enum):
    HBOND_LOSS = "hbond_loss"
    MG_BRIDGE_DISRUPTION = "mg_bridge_disruption"
    CHARGE_CHANGE = "charge_change"
    OTHER = "other"


def _enum(py_enum: type[enum.Enum], name: str) -> SAEnum:
    return SAEnum(
        py_enum,
        name=name,
        native_enum=False,
        create_constraint=True,
        validate_strings=True,
        values_callable=lambda e: [m.value for m in e],
    )


# --------------------------------------------------------------------- lookup tables
class Species(Base):
    __tablename__ = "species"

    species_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    scientific_name: Mapped[str] = mapped_column(String(120), unique=True, nullable=False)
    ncbi_taxid: Mapped[Optional[int]] = mapped_column(Integer)
    short_code: Mapped[Optional[str]] = mapped_column(String(10))

    isolates: Mapped[List["Isolate"]] = relationship(back_populates="species")
    positions: Mapped[List["TargetPosition"]] = relationship(back_populates="species")


class Habitat(Base):
    __tablename__ = "habitats"

    habitat_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(60), unique=True, nullable=False)

    isolates: Mapped[List["Isolate"]] = relationship(back_populates="habitat")


class AminoAcidProperty(Base):
    """Physicochemical lookup. Hydropathy = Eisenberg consensus scale (1984)."""

    __tablename__ = "amino_acid_properties"

    aa: Mapped[str] = mapped_column(String(1), primary_key=True)
    name: Mapped[str] = mapped_column(String(20), nullable=False)
    hydropathy_value: Mapped[float] = mapped_column(Float, nullable=False)
    charge: Mapped[int] = mapped_column(SmallInteger, nullable=False, default=0)
    volume: Mapped[Optional[float]] = mapped_column(Float)  # cubic angstrom
    polarity: Mapped[Optional[str]] = mapped_column(String(20))
    has_hydroxyl: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    scale_name: Mapped[str] = mapped_column(
        String(80), nullable=False, default="Eisenberg consensus (1984)"
    )


class TargetPosition(Base):
    """Species-specific numbering of the five codon positions."""

    __tablename__ = "target_positions"
    __table_args__ = (UniqueConstraint("species_id", "gene", "residue_number"),)

    position_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    gene: Mapped[str] = mapped_column(String(10), nullable=False)  # gyrA, parC, recN
    residue_number: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    label: Mapped[str] = mapped_column(String(20), nullable=False)  # e.g. GyrA-83
    species_id: Mapped[int] = mapped_column(ForeignKey("species.species_id"), nullable=False)
    reference_accession: Mapped[Optional[str]] = mapped_column(String(40))

    species: Mapped[Species] = relationship(back_populates="positions")


# ------------------------------------------------------------------------- septatypes
class Septatype(Base):
    __tablename__ = "septatypes"

    septatype_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    septatype_string: Mapped[str] = mapped_column(String(30), unique=True, nullable=False)

    gyra83_aa: Mapped[str] = mapped_column(String(1), nullable=False)
    gyra87_aa: Mapped[str] = mapped_column(String(1), nullable=False)
    parc80_aa: Mapped[str] = mapped_column(String(1), nullable=False)
    parc84_aa: Mapped[str] = mapped_column(String(1), nullable=False)
    recn346_aa: Mapped[str] = mapped_column(String(1), nullable=False)
    qnrb2: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    qnrs1: Mapped[int] = mapped_column(SmallInteger, nullable=False)

    mean_hydrophobicity: Mapped[Optional[float]] = mapped_column(Float)  # <H> GyrA-ParC
    predicted_phenotype: Mapped[Optional[SeptaPrediction]] = mapped_column(
        _enum(SeptaPrediction, "septa_prediction")
    )
    mechanism_class: Mapped[Optional[str]] = mapped_column(String(80))
    is_novel: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    first_seen_date: Mapped[Optional[dt.date]] = mapped_column(Date)
    isolate_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    isolates: Mapped[List["Isolate"]] = relationship(back_populates="septatype")
    annotations: Mapped[List["MechanismAnnotation"]] = relationship(back_populates="septatype")

    __table_args__ = (
        CheckConstraint("qnrb2 IN (0,1) AND qnrs1 IN (0,1)", name="qnr_binary"),
    )


class MechanismAnnotation(Base):
    """Feeds the Module D interactive structure viewer."""

    __tablename__ = "mechanism_annotations"

    annotation_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    septatype_id: Mapped[Optional[int]] = mapped_column(ForeignKey("septatypes.septatype_id"))
    residue_combination: Mapped[Optional[str]] = mapped_column(String(60))
    effect_type: Mapped[EffectType] = mapped_column(_enum(EffectType, "effect_type"), nullable=False)
    description: Mapped[str] = mapped_column(Text, nullable=False)
    pdb_id: Mapped[Optional[str]] = mapped_column(String(10))
    residue_pairs: Mapped[Optional[dict]] = mapped_column(JSONType)
    viewer_state: Mapped[Optional[dict]] = mapped_column(JSONType)

    septatype: Mapped[Optional[Septatype]] = relationship(back_populates="annotations")


# ----------------------------------------------------------------- users / submissions
class User(Base):
    __tablename__ = "users"

    user_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    email: Mapped[str] = mapped_column(String(255), unique=True, nullable=False)
    hashed_password: Mapped[Optional[str]] = mapped_column(String(255))
    role: Mapped[UserRole] = mapped_column(
        _enum(UserRole, "user_role"), nullable=False, default=UserRole.SUBMITTER
    )
    created_at: Mapped[dt.datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), nullable=False
    )


class AnalysisJob(Base):
    __tablename__ = "analysis_jobs"

    job_id: Mapped[str] = mapped_column(String(36), primary_key=True)  # UUID4 string
    input_type: Mapped[JobInputType] = mapped_column(_enum(JobInputType, "job_input_type"), nullable=False)
    status: Mapped[JobStatus] = mapped_column(
        _enum(JobStatus, "job_status"), nullable=False, default=JobStatus.QUEUED
    )
    result: Mapped[Optional[dict]] = mapped_column(JSONType)
    error: Mapped[Optional[str]] = mapped_column(Text)
    created_at: Mapped[dt.datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), nullable=False
    )
    expires_at: Mapped[Optional[dt.datetime]] = mapped_column(DateTime(timezone=True))


class Submission(Base):
    __tablename__ = "submissions"

    submission_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    submitter_id: Mapped[Optional[int]] = mapped_column(ForeignKey("users.user_id"))
    status: Mapped[SubmissionStatus] = mapped_column(
        _enum(SubmissionStatus, "submission_status"),
        nullable=False,
        default=SubmissionStatus.PENDING,
        index=True,
    )

    # submitter-supplied metadata
    country: Mapped[Optional[str]] = mapped_column(String(80))
    collection_year: Mapped[Optional[int]] = mapped_column(SmallInteger)
    source: Mapped[Optional[str]] = mapped_column(String(80))
    experimental_mic: Mapped[Optional[float]] = mapped_column(Numeric(10, 4))
    mic_unit: Mapped[Optional[str]] = mapped_column(String(20))
    experimental_phenotype: Mapped[Optional[Phenotype]] = mapped_column(_enum(Phenotype, "phenotype"))

    # who and what was submitted (Phase 5)
    contact_email: Mapped[Optional[str]] = mapped_column(String(255))
    submitter_name: Mapped[Optional[str]] = mapped_column(String(120))
    isolate_name: Mapped[Optional[str]] = mapped_column(String(200))
    declared_species: Mapped[Optional[str]] = mapped_column(String(60))
    habitat: Mapped[Optional[str]] = mapped_column(String(60))
    country_iso3: Mapped[Optional[str]] = mapped_column(String(3))
    mic_method: Mapped[Optional[str]] = mapped_column(String(60))
    file_sha256: Mapped[Optional[str]] = mapped_column(String(64), index=True)
    status_token_hash: Mapped[Optional[str]] = mapped_column(String(64))

    # BV-BRC lookup
    bv_brc_id_supplied: Mapped[Optional[str]] = mapped_column(String(30), index=True)
    bv_brc_name_fetched: Mapped[Optional[str]] = mapped_column(Text)
    bv_brc_lookup_status: Mapped[Optional[str]] = mapped_column(String(15))  # not_provided | ok | not_found | error
    bvbrc_amr: Mapped[Optional[dict]] = mapped_column(JSONType)  # {prediction, raw, method} for benchmarking

    # background analysis (Phase 5)
    analysis_state: Mapped[Optional[str]] = mapped_column(String(12), index=True)  # queued | running | done | failed
    analysis_error: Mapped[Optional[str]] = mapped_column(Text)
    prediction_result: Mapped[Optional[dict]] = mapped_column(JSONType)  # full engine output
    predicted_septatype_string: Mapped[Optional[str]] = mapped_column(String(30))
    septatyper_prediction: Mapped[Optional[str]] = mapped_column(String(15))  # Cres | Csus | Indeterminate
    flags: Mapped[Optional[list]] = mapped_column(JSONType)  # list of flag codes

    # analysis
    file_path: Mapped[Optional[str]] = mapped_column(Text)
    job_id: Mapped[Optional[str]] = mapped_column(ForeignKey("analysis_jobs.job_id"))
    predicted_septatype_id: Mapped[Optional[int]] = mapped_column(ForeignKey("septatypes.septatype_id"))
    is_novel_septatype: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    prediction_matches_phenotype: Mapped[Optional[bool]] = mapped_column(Boolean)

    # curation
    curator_id: Mapped[Optional[int]] = mapped_column(ForeignKey("users.user_id"))
    curator_notes: Mapped[Optional[str]] = mapped_column(Text)
    decided_at: Mapped[Optional[dt.datetime]] = mapped_column(DateTime(timezone=True))
    created_at: Mapped[dt.datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), nullable=False
    )

    isolate: Mapped[Optional["Isolate"]] = relationship(back_populates="submission", uselist=False)


# -------------------------------------------------------------------------- isolates
_SEPTA_CORRECT = (
    "((septatyper_prediction = 'Cres' AND lab_phenotype = 'Resistant') "
    "OR (septatyper_prediction = 'Csus' AND lab_phenotype = 'Susceptible'))"
)


class Isolate(Base):
    __tablename__ = "isolates"

    bv_brc_genome_id: Mapped[str] = mapped_column(String(30), primary_key=True)
    official_genome_name: Mapped[Optional[str]] = mapped_column(Text)
    isolate_name: Mapped[Optional[str]] = mapped_column(String(200))

    species_id: Mapped[int] = mapped_column(ForeignKey("species.species_id"), nullable=False, index=True)
    habitat_id: Mapped[Optional[int]] = mapped_column(ForeignKey("habitats.habitat_id"), index=True)
    mlst_st: Mapped[Optional[str]] = mapped_column(String(30), index=True)

    collection_year_start: Mapped[Optional[int]] = mapped_column(SmallInteger)
    collection_year_end: Mapped[Optional[int]] = mapped_column(SmallInteger)
    collection_year_raw: Mapped[Optional[str]] = mapped_column(String(20))
    country_name: Mapped[Optional[str]] = mapped_column(String(80))
    country_iso3: Mapped[Optional[str]] = mapped_column(String(3), index=True)

    lab_phenotype: Mapped[Phenotype] = mapped_column(_enum(Phenotype, "phenotype"), nullable=False, index=True)
    mic_value: Mapped[Optional[float]] = mapped_column(Numeric(10, 4))
    mic_unit: Mapped[Optional[str]] = mapped_column(String(20))
    mic_method: Mapped[Optional[str]] = mapped_column(String(60))

    # BV-BRC computational benchmark (raw value and method kept for provenance)
    bv_brc_computational_prediction: Mapped[Optional[BvbrcPrediction]] = mapped_column(
        _enum(BvbrcPrediction, "bvbrc_prediction")
    )
    bv_brc_prediction_raw: Mapped[Optional[str]] = mapped_column(Text)
    bv_brc_prediction_method: Mapped[Optional[str]] = mapped_column(Text)
    bv_brc_fetched_at: Mapped[Optional[dt.datetime]] = mapped_column(DateTime(timezone=True))

    septatyper_prediction: Mapped[Optional[SeptaPrediction]] = mapped_column(
        _enum(SeptaPrediction, "septa_prediction")
    )
    septatype_id: Mapped[Optional[int]] = mapped_column(ForeignKey("septatypes.septatype_id"), index=True)

    # Generated benchmarking flags (NULL until the relevant prediction exists).
    # PostgreSQL forbids generated columns referencing other generated columns,
    # so each expression is written out in full.
    septatyper_correct: Mapped[Optional[bool]] = mapped_column(
        Boolean,
        Computed(
            f"CASE WHEN septatyper_prediction IS NULL THEN NULL "
            f"WHEN {_SEPTA_CORRECT} THEN TRUE ELSE FALSE END",
            persisted=True,
        ),
    )
    bvbrc_correct: Mapped[Optional[bool]] = mapped_column(
        Boolean,
        Computed(
            "CASE WHEN bv_brc_computational_prediction IS NULL THEN NULL "
            "WHEN bv_brc_computational_prediction = lab_phenotype THEN TRUE ELSE FALSE END",
            persisted=True,
        ),
    )
    septatyper_wins: Mapped[Optional[bool]] = mapped_column(
        Boolean,
        Computed(
            "CASE WHEN septatyper_prediction IS NULL "
            "OR bv_brc_computational_prediction IS NULL THEN NULL "
            f"WHEN {_SEPTA_CORRECT} "
            "AND bv_brc_computational_prediction <> lab_phenotype THEN TRUE ELSE FALSE END",
            persisted=True,
        ),
        index=True,
    )

    source: Mapped[IsolateSource] = mapped_column(
        _enum(IsolateSource, "isolate_source"), nullable=False, default=IsolateSource.FOUNDING
    )
    submission_id: Mapped[Optional[int]] = mapped_column(ForeignKey("submissions.submission_id"))
    date_added: Mapped[dt.datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), nullable=False
    )

    species: Mapped[Species] = relationship(back_populates="isolates")
    habitat: Mapped[Optional[Habitat]] = relationship(back_populates="isolates")
    septatype: Mapped[Optional[Septatype]] = relationship(back_populates="isolates")
    submission: Mapped[Optional[Submission]] = relationship(back_populates="isolate")
    genotype: Mapped[Optional["IsolateGenotype"]] = relationship(
        back_populates="isolate", uselist=False, cascade="all, delete-orphan"
    )
    pmqr_calls: Mapped[List["PMQRCall"]] = relationship(
        back_populates="isolate", cascade="all, delete-orphan"
    )

    __table_args__ = (
        Index("ix_isolates_year_range", "collection_year_start", "collection_year_end"),
    )


class IsolateGenotype(Base):
    """1:1 with isolates. Wide layout keeps septatype filtering fast.

    Codon columns are filled by the BLASTn pipeline (Phase 3). Table S1 only
    provides amino acids, so *_codon stays NULL after the initial ingestion.
    """

    __tablename__ = "isolate_genotypes"

    bv_brc_genome_id: Mapped[str] = mapped_column(
        ForeignKey("isolates.bv_brc_genome_id", ondelete="CASCADE"), primary_key=True
    )

    gyra83_codon: Mapped[Optional[str]] = mapped_column(String(3))
    gyra83_aa: Mapped[Optional[str]] = mapped_column(String(1))
    gyra87_codon: Mapped[Optional[str]] = mapped_column(String(3))
    gyra87_aa: Mapped[Optional[str]] = mapped_column(String(1))
    parc80_codon: Mapped[Optional[str]] = mapped_column(String(3))
    parc80_aa: Mapped[Optional[str]] = mapped_column(String(1))
    parc84_codon: Mapped[Optional[str]] = mapped_column(String(3))
    parc84_aa: Mapped[Optional[str]] = mapped_column(String(1))
    recn346_codon: Mapped[Optional[str]] = mapped_column(String(3))
    recn346_aa: Mapped[Optional[str]] = mapped_column(String(1))
    recn346_synonymous_flag: Mapped[Optional[bool]] = mapped_column(Boolean)  # CTA vs CTG in S. enterica

    alignment_identity: Mapped[Optional[float]] = mapped_column(Float)
    coverage: Mapped[Optional[float]] = mapped_column(Float)
    reference_used: Mapped[Optional[str]] = mapped_column(String(80))
    pipeline_version: Mapped[Optional[str]] = mapped_column(String(40))

    isolate: Mapped[Isolate] = relationship(back_populates="genotype")


class PMQRCall(Base):
    __tablename__ = "pmqr_calls"
    __table_args__ = (UniqueConstraint("bv_brc_genome_id", "gene"),)

    call_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    bv_brc_genome_id: Mapped[str] = mapped_column(
        ForeignKey("isolates.bv_brc_genome_id", ondelete="CASCADE"), nullable=False, index=True
    )
    gene: Mapped[str] = mapped_column(String(20), nullable=False)  # qnrB2, qnrS1, ...
    present: Mapped[bool] = mapped_column(Boolean, nullable=False)
    identity: Mapped[Optional[float]] = mapped_column(Float)
    coverage: Mapped[Optional[float]] = mapped_column(Float)
    contig_id: Mapped[Optional[str]] = mapped_column(String(120))

    isolate: Mapped[Isolate] = relationship(back_populates="pmqr_calls")


class AuditLog(Base):
    __tablename__ = "audit_log"

    log_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    user_id: Mapped[Optional[int]] = mapped_column(ForeignKey("users.user_id"))
    action: Mapped[str] = mapped_column(String(60), nullable=False)
    target_table: Mapped[Optional[str]] = mapped_column(String(60))
    target_id: Mapped[Optional[str]] = mapped_column(String(60))
    detail: Mapped[Optional[dict]] = mapped_column(JSONType)
    timestamp: Mapped[dt.datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), nullable=False
    )
