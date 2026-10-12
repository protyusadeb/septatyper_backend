"""Idempotent seeding of lookup tables."""
from __future__ import annotations

from sqlalchemy import select
from sqlalchemy.orm import Session

from .models import AminoAcidProperty, Habitat, Species, TargetPosition

# (name, ncbi_taxid, short_code)
SPECIES = [
    ("Escherichia coli", 562, "Eco"),
    ("Shigella sonnei", 624, "Sso"),
    ("Salmonella enterica subspecies I", 28901, "Sen"),
    ("Klebsiella pneumoniae", 573, "Kpn"),
]

HABITATS = ["Fecal", "Urine", "Respiratory"]

# aa: (name, Eisenberg consensus hydrophobicity, charge, volume A^3, polarity class, has -OH)
# Eisenberg, Schwarz, Komaromy & Wall 1984, J Mol Biol 179:125-142 (original consensus scale; the
# manuscript quotes Ser -0.18, Ile 1.38, Leu 1.06, Phe 1.19). The 20 values sum to zero, so a mean
# of 0 is the scale average. Volumes: Zamyatnin (1972).
SCALE_NAME = "Eisenberg consensus (1984)"
# Polarity classes are a convention for browsing only; the Phase 3 rule engine
# applies the grammar's own definitions.
AMINO_ACIDS = {
    "A": ("Alanine", 0.62, 0, 88.6, "hydrophobic", False),
    "R": ("Arginine", -2.53, 1, 173.4, "positive", False),
    "N": ("Asparagine", -0.78, 0, 114.1, "polar", False),
    "D": ("Aspartate", -0.9, -1, 111.1, "negative", False),
    "C": ("Cysteine", 0.29, 0, 108.5, "hydrophobic", False),
    "Q": ("Glutamine", -0.85, 0, 143.8, "polar", False),
    "E": ("Glutamate", -0.74, -1, 138.4, "negative", False),
    "G": ("Glycine", 0.48, 0, 60.1, "hydrophobic", False),
    "H": ("Histidine", -0.4, 0, 153.2, "positive", False),
    "I": ("Isoleucine", 1.38, 0, 166.7, "hydrophobic", False),
    "L": ("Leucine", 1.06, 0, 166.7, "hydrophobic", False),
    "K": ("Lysine", -1.5, 1, 168.6, "positive", False),
    "M": ("Methionine", 0.64, 0, 162.9, "hydrophobic", False),
    "F": ("Phenylalanine", 1.19, 0, 189.9, "hydrophobic", False),
    "P": ("Proline", 0.12, 0, 112.7, "hydrophobic", False),
    "S": ("Serine", -0.18, 0, 89.0, "polar", True),
    "T": ("Threonine", -0.05, 0, 116.1, "polar", True),
    "W": ("Tryptophan", 0.81, 0, 227.8, "hydrophobic", False),
    "Y": ("Tyrosine", 0.26, 0, 193.6, "polar", True),
    "V": ("Valine", 1.08, 0, 140.0, "hydrophobic", False),
}

# (gene, residue number, label). Numbering follows the Table S1 column labels.
# Species-specific offsets (notably RecN-346 in Salmonella) must be confirmed
# against the reference sequences in Phase 3.
POSITIONS = [
    ("gyrA", 83, "GyrA-83"),
    ("gyrA", 87, "GyrA-87"),
    ("parC", 80, "ParC-80"),
    ("parC", 84, "ParC-84"),
    ("recN", 346, "RecN-346"),
]


def seed_lookups(session: Session) -> None:
    for name, taxid, code in SPECIES:
        if not session.scalar(select(Species).where(Species.scientific_name == name)):
            session.add(Species(scientific_name=name, ncbi_taxid=taxid, short_code=code))
    for name in HABITATS:
        if not session.scalar(select(Habitat).where(Habitat.name == name)):
            session.add(Habitat(name=name))
    for aa, (name, hyd, charge, vol, pol, oh) in AMINO_ACIDS.items():
        row = session.get(AminoAcidProperty, aa)
        if row is None:
            session.add(
                AminoAcidProperty(
                    aa=aa, name=name, hydropathy_value=hyd, charge=charge,
                    volume=vol, polarity=pol, has_hydroxyl=oh, scale_name=SCALE_NAME,
                )
            )
        else:  # re-running the seed updates values seeded with an earlier scale
            row.hydropathy_value, row.scale_name = hyd, SCALE_NAME
    session.flush()
    for sp in session.scalars(select(Species)).all():
        for gene, num, label in POSITIONS:
            exists = session.scalar(
                select(TargetPosition).where(
                    TargetPosition.species_id == sp.species_id,
                    TargetPosition.gene == gene,
                    TargetPosition.residue_number == num,
                )
            )
            if not exists:
                session.add(
                    TargetPosition(gene=gene, residue_number=num, label=label, species_id=sp.species_id)
                )
    session.commit()
