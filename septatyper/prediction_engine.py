"""SeptaTyper prediction engine (Module A).

Pipeline for one uploaded genome assembly (nucleotide FASTA):

1. Sanitise the FASTA and build a temporary BLAST nucleotide database from it.
2. BLASTn the local reference genes (gyrA, parC, recN) against the genome.
3. Pull the five target codons out of the alignments at fixed, gene-relative
   coordinates (codon first, translation second).
4. BLASTn the PMQR references (qnrB2, qnrS1); present only if identity >= 80%
   AND coverage >= 90%.
5. Run the two-tier grammar: Tier 1 physicochemical baseline, Tier 2 overrides.

Reference layout (see references/README.md):
    references/gyrA.fasta  parC.fasta  recN.fasta     CDS sequences starting at the ATG
    references/pmqr/qnrB2.fasta  qnrS1.fasta           any number of alleles per file
Headers of the chromosomal files: ">gene|Species name|accession". The species field
drives species inference.

Command line:
    python -m septatyper.prediction_engine --check-refs
    python -m septatyper.prediction_engine genome.fasta [--species "Shigella sonnei"]
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from itertools import product
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .seed import AMINO_ACIDS

log = logging.getLogger("septatyper.engine")

PIPELINE_VERSION = "septatyper-engine/1.0"

# --------------------------------------------------------------------------- constants
ECOLI = "Escherichia coli"
SSONNEI = "Shigella sonnei"
SENTERICA = "Salmonella enterica subspecies I"
KPNEU = "Klebsiella pneumoniae"
SPECIES = (ECOLI, SSONNEI, SENTERICA, KPNEU)

# Gene-relative, 1-based coordinate of the FIRST nucleotide of each target codon.
# These come from the validated Perl workflow and equal 3 * codon_number - 2.
TARGETS: Dict[str, Tuple[str, int, int]] = {
    # key: (gene, codon number, first nucleotide)
    "gyrA_83": ("gyrA", 83, 247),   # nt 247-249
    "gyrA_87": ("gyrA", 87, 259),   # nt 259-261
    "parC_80": ("parC", 80, 238),   # nt 238-240
    "parC_84": ("parC", 84, 250),   # nt 250-252
    "recN_346": ("recN", 346, 1036),  # nt 1036-1038
}
assert all(start == 3 * codon - 2 for _, codon, start in TARGETS.values())

TARGET_LABELS = {
    "gyrA_83": "GyrA-83", "gyrA_87": "GyrA-87", "parC_80": "ParC-80",
    "parC_84": "ParC-84", "recN_346": "RecN-346",
}
# Residues expected in a wild-type reference; used only for reference sanity checks.
EXPECTED_WT = {"gyrA_83": "S", "gyrA_87": "D", "parC_80": "S", "parC_84": "E", "recN_346": "L"}
GENES = ("gyrA", "parC", "recN")
PMQR_GENES = ("qnrB2", "qnrS1")

PMQR_MIN_IDENTITY = 80.0   # percent nucleotide identity
PMQR_MIN_COVERAGE = 90.0   # percent of the reference covered
MIN_GENE_IDENTITY = 85.0   # chromosomal gene hits below this are ignored

BLAST_FIELDS = (
    "qseqid sseqid pident length nident qlen qstart qend sstart send evalue bitscore qseq sseq"
)

EISENBERG: Dict[str, float] = {aa: v[1] for aa, v in AMINO_ACIDS.items()}
HYDROPHOBICITY_SCALE = "Eisenberg consensus (1984)"

# Class definitions used by the grammar. Tyrosine counts as polar (-OH), as specified.
HYDROPHOBIC = frozenset("ACFGILMPVW")
POLAR_OH = frozenset("STY")
NEGATIVE = frozenset("DE")
CHARGED = frozenset("DEKRH")
H_RESCUE_THRESHOLD = 0.0  # <H> (mean of the 4 residues) at or above the scale mean (0) can tip a single-vote tetratype
# S. sonnei tetratypes with a single resistance vote that are nonetheless Resistant (Tier 2).
SONNEI_SINGLE_VOTE_TETRATYPES = frozenset({"L-D-S-E", "S-Y-S-E"})
# Manuscript notation: RecN-346 is written Lcta / Lctg (instead of plain L) only for these
# species (S. enterica and K. pneumoniae) when the tetratype is F-D-S-E. Extend the set if the notation should widen.
CODON_LABELLED_SPECIES = frozenset({SENTERICA, KPNEU})
CODON_LABELLED_TETRATYPE = "F-D-S-E"

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(product(_BASES, repeat=3))}

IUPAC_NT = set("ACGTRYSWKMBDHVN")


class EngineError(RuntimeError):
    """Configuration or runtime failure (missing BLAST, bad references, BLAST crash)."""


class InputError(ValueError):
    """The uploaded sequence is unusable."""


# ----------------------------------------------------------------------- species helper
def canonical_species(label: Optional[str]) -> Optional[str]:
    """Map free-text species names or short codes to the four canonical names."""
    if not label:
        return None
    s = re.sub(r"[_\s]+", " ", label.strip().lower())
    table = {
        "escherichia coli": ECOLI, "e coli": ECOLI, "eco": ECOLI,
        "shigella sonnei": SSONNEI, "s sonnei": SSONNEI, "sso": SSONNEI,
        "klebsiella pneumoniae": KPNEU, "k pneumoniae": KPNEU, "kpn": KPNEU,
        "sen": SENTERICA, "s enterica": SENTERICA,
    }
    if s in table:
        return table[s]
    if s.startswith("salmonella enterica"):
        return SENTERICA
    return None


# ------------------------------------------------------------------------------ FASTA
def read_fasta(text: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, List[str]]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            records.append((line[1:].strip(), []))
        elif records:
            records[-1][1].append(line)
        else:
            raise InputError("FASTA must start with a '>' header line")
    return [(h, "".join(parts)) for h, parts in records]


def translate(codon: str) -> Optional[str]:
    return CODON_TABLE.get(codon.upper())


# ---------------------------------------------------------------------------- BLAST io
@dataclass
class HSP:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    nident: int
    qlen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    bitscore: float
    qseq: str
    sseq: str

    @property
    def gene(self) -> str:
        return self.qseqid.split("|", 1)[0]

    @property
    def ref_idx(self) -> int:
        return int(self.qseqid.split("|", 1)[1])

    @property
    def strand(self) -> str:
        return "+" if self.sstart <= self.send else "-"


def parse_blast_tabular(text: str) -> List[HSP]:
    hsps: List[HSP] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) != 14:
            raise EngineError(f"Unexpected BLAST output line with {len(f)} fields")
        hsps.append(
            HSP(
                qseqid=f[0], sseqid=f[1], pident=float(f[2]), length=int(f[3]), nident=int(f[4]),
                qlen=int(f[5]), qstart=int(f[6]), qend=int(f[7]), sstart=int(f[8]), send=int(f[9]),
                bitscore=float(f[11]), qseq=f[12].upper(), sseq=f[13].upper(),
            )
        )
    return hsps


def extract_codon(hsp: HSP, ref_start: int) -> Tuple[Optional[str], str]:
    """Return (codon, status) for reference positions ref_start..ref_start+2.

    The query is the REFERENCE gene, so qstart/qend are reference coordinates and
    sseq holds the genome bases in the reference orientation (BLAST reverse-complements
    minus-strand hits in the alignment), so no strand handling is needed here.
    """
    end = ref_start + 2
    if hsp.qstart > ref_start or hsp.qend < end:
        return None, "not_covered"
    pos = hsp.qstart
    bases: List[str] = []
    for qc, sc in zip(hsp.qseq, hsp.sseq):
        if qc == "-":  # insertion in the genome relative to the reference
            if bases and len(bases) < 3:
                return None, "indel_in_codon"
            continue
        if ref_start <= pos <= end:
            if sc == "-":  # deletion in the genome
                return None, "indel_in_codon"
            bases.append(sc)
        pos += 1
        if pos > end:
            break
    if len(bases) != 3:
        return None, "not_covered"
    codon = "".join(bases)
    if any(b not in "ACGT" for b in codon):
        return None, "ambiguous_base"
    return codon, "ok"


def _union_length(intervals: Iterable[Tuple[int, int]]) -> int:
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1  # type: ignore[operator]
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1  # type: ignore[operator]
    return total


def summarise_hsps(hsps: Sequence[HSP]) -> Tuple[float, float]:
    """(identity %, coverage % of the reference) for a set of HSPs of one reference."""
    qlen = hsps[0].qlen
    cov = 100.0 * _union_length((h.qstart, h.qend) for h in hsps) / qlen
    ident = 100.0 * sum(h.nident for h in hsps) / max(1, sum(h.length for h in hsps))
    return round(ident, 2), round(min(cov, 100.0), 2)


# --------------------------------------------------------------------------- grammar
@dataclass
class Tier1Result:
    baseline: str                     # "Resistant" | "Susceptible"
    votes: int
    decided_by: str                   # "class_votes" | "hydrophobicity"
    mean_hydrophobicity: float
    detail: List[dict]


@dataclass
class Tier2Result:
    final: str
    override_applied: bool
    rules: List[str]


def mean_hydrophobicity(aa: Mapping[str, Optional[str]]) -> Optional[float]:
    """<H> GyrA-ParC: mean Eisenberg value of GyrA-83, GyrA-87, ParC-80, ParC-84."""
    keys = ("gyrA_83", "gyrA_87", "parC_80", "parC_84")
    residues = [aa.get(k) for k in keys]
    if any(r is None or r not in EISENBERG for r in residues):
        return None
    return sum(EISENBERG[r] for r in residues) / 4.0  # type: ignore[index]


def tier1(aa: Mapping[str, Optional[str]]) -> Tier1Result:
    """Tier 1: physicochemical baseline from the four quinolone-target residues.

    Each position casts a resistance vote:
      GyrA-83  hydrophobic = 1 ; polar with -OH (S, T, Y) = 0
      GyrA-87  loss of negative charge (anything but D/E) = 1 ; D/E = 0
      ParC-80  hydrophobic = 1 ; polar with -OH = 0
      ParC-84  uncharged = 1 ; K/R/H = 0 ; D/E = 1 ONLY when ParC-80 is hydrophobic
               (cooperative: S-E keeps the water-metal ion bridge, I-E disrupts it)
    Two or more votes = Resistant, none = Susceptible. A single vote is resolved by
    <H>: at or above H_RESCUE_THRESHOLD (the scale mean) the hydrophobic
    environment tips it to Resistant.
    """
    g83, g87, p80, p84 = (aa["gyrA_83"], aa["gyrA_87"], aa["parC_80"], aa["parC_84"])
    detail: List[dict] = []

    def add(pos: str, res: str, cls: str, vote: int) -> int:
        detail.append({"position": pos, "residue": res, "class": cls, "resistance_vote": vote})
        return vote

    v = 0
    v += add("GyrA-83", g83, "hydrophobic" if g83 in HYDROPHOBIC else "polar_OH" if g83 in POLAR_OH else "other",
             int(g83 in HYDROPHOBIC))
    v += add("GyrA-87", g87, "negative" if g87 in NEGATIVE else "polar_or_uncharged",
             int(g87 not in NEGATIVE))
    v += add("ParC-80", p80, "hydrophobic" if p80 in HYDROPHOBIC else "polar_OH" if p80 in POLAR_OH else "other",
             int(p80 in HYDROPHOBIC))
    if p84 in NEGATIVE:
        v += add("ParC-84", p84, "negative_with_hydrophobic_ParC80" if p80 in HYDROPHOBIC else "negative",
                 int(p80 in HYDROPHOBIC))
    elif p84 in CHARGED:
        v += add("ParC-84", p84, "positive_or_charged", 0)
    else:
        v += add("ParC-84", p84, "uncharged", 1)

    h = mean_hydrophobicity(aa)
    assert h is not None
    if v >= 2:
        baseline, by = "Resistant", "class_votes"
    elif v == 0:
        baseline, by = "Susceptible", "class_votes"
    elif h >= H_RESCUE_THRESHOLD:
        baseline, by = "Resistant", "hydrophobicity"
    else:
        baseline, by = "Susceptible", "hydrophobicity"
    return Tier1Result(baseline, v, by, h, detail)


def tier2(
    baseline: str,
    aa: Mapping[str, Optional[str]],
    pmqr: Mapping[str, int],
    species: Optional[str],
    recn_codon: Optional[str],
) -> Tier2Result:
    """Tier 2 overrides. They only ever upgrade a Susceptible baseline to Resistant.

    qnr_present                  qnrB2 or qnrS1 present (any species)
    sonnei_recN_L346F            S. sonnei, S-D-S-E, RecN-346 F
    enterica_recN_ctg            S. enterica, F-D-S-E, RecN-346 codon CTG
    sonnei_single_vote_override  S. sonnei, L-D-S-E or S-Y-S-E
    """
    if baseline == "Resistant":
        return Tier2Result("Resistant", False, [])
    tetratype = "-".join(str(aa[k]) for k in ("gyrA_83", "gyrA_87", "parC_80", "parC_84"))
    rules: List[str] = []
    if pmqr.get("qnrB2") == 1 or pmqr.get("qnrS1") == 1:
        rules.append("qnr_present")
    if species == SSONNEI and tetratype == "S-D-S-E" and aa.get("recN_346") == "F":
        rules.append("sonnei_recN_L346F")
    if species == SENTERICA and tetratype == "F-D-S-E" and recn_codon == "CTG":
        rules.append("enterica_recN_ctg")
    if species == SSONNEI and tetratype in SONNEI_SINGLE_VOTE_TETRATYPES:
        rules.append("sonnei_single_vote_override")
    return Tier2Result("Resistant" if rules else "Susceptible", bool(rules), rules)


def aa_level_septatype(septatype: Optional[str]) -> Optional[str]:
    """Strip the codon notation: F-D-S-E-Lctg-0-1 -> F-D-S-E-L-0-1."""
    return re.sub(r"-L(?:cta|ctg)-", "-L-", septatype) if septatype else None


def recn_label(aa: Mapping[str, Optional[str]], species: Optional[str], recn_codon: Optional[str]) -> Optional[str]:
    """RecN-346 symbol for the septatype string: Lcta/Lctg for S. enterica F-D-S-E, else the amino acid."""
    residue = aa.get("recN_346")
    tetratype = "-".join(str(aa.get(k)) for k in ("gyrA_83", "gyrA_87", "parC_80", "parC_84"))
    if (species in CODON_LABELLED_SPECIES and tetratype == CODON_LABELLED_TETRATYPE
            and residue == "L" and recn_codon in ("CTA", "CTG")):
        return "L" + recn_codon.lower()
    return residue


@dataclass
class GrammarResult:
    septatype: Optional[str]
    prediction: str                    # "Cres" | "Csus" | "Indeterminate"
    phenotype: Optional[str]
    mean_hydrophobicity: Optional[float]
    tier1: Optional[Tier1Result]
    tier2: Optional[Tier2Result]
    septatype_aa_level: Optional[str] = None


def run_grammar(
    aa: Mapping[str, Optional[str]],
    pmqr: Mapping[str, int],
    species: Optional[str],
    recn_codon: Optional[str],
) -> GrammarResult:
    five = [aa.get(k) for k in TARGETS]
    septatype = None
    if all(five):
        parts = list(five[:4]) + [recn_label(aa, species, recn_codon)] + [str(pmqr["qnrB2"]), str(pmqr["qnrS1"])]
        septatype = "-".join(parts)  # type: ignore[arg-type]
    plain = aa_level_septatype(septatype)
    if any(aa.get(k) not in EISENBERG for k in ("gyrA_83", "gyrA_87", "parC_80", "parC_84")):
        return GrammarResult(septatype, "Indeterminate", None, None, None, None, plain)
    t1 = tier1(aa)
    t2 = tier2(t1.baseline, aa, pmqr, species, recn_codon)
    return GrammarResult(
        septatype, "Cres" if t2.final == "Resistant" else "Csus", t2.final, t1.mean_hydrophobicity, t1, t2, plain
    )


# ------------------------------------------------------------------------- references
@dataclass
class Ref:
    gene: str
    label: str
    species: Optional[str]    # canonical, or None when unknown
    species_raw: str
    seq: str


@dataclass
class ReferenceReport:
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def _load_refs(ref_dir: Path) -> Tuple[Dict[str, List[Ref]], ReferenceReport]:
    rep = ReferenceReport()
    refs: Dict[str, List[Ref]] = {}
    for gene in GENES + PMQR_GENES:
        path = ref_dir / f"{gene}.fasta" if gene in GENES else ref_dir / "pmqr" / f"{gene}.fasta"
        if not path.exists():
            rep.errors.append(f"Missing reference file: {path}")
            continue
        recs = read_fasta(path.read_text())
        if not recs:
            rep.errors.append(f"No sequences in {path}")
        out = []
        for header, seq in recs:
            parts = header.split("|")
            raw = parts[1].strip() if len(parts) > 1 else ""
            out.append(Ref(gene, header, canonical_species(raw), raw, seq.upper()))
        refs[gene] = out
    return refs, rep


def validate_references(refs: Mapping[str, List[Ref]], rep: Optional[ReferenceReport] = None) -> ReferenceReport:
    """Check the assumptions the fixed coordinates depend on."""
    rep = rep or ReferenceReport()
    for key, (gene, codon_no, start) in TARGETS.items():
        for r in refs.get(gene, []):
            tag = f"{r.label}"
            if not r.seq.startswith("ATG"):
                rep.errors.append(f"{tag}: does not start with ATG; fixed coordinates would be shifted")
                continue
            if len(r.seq) % 3 != 0:
                rep.errors.append(f"{tag}: length {len(r.seq)} is not a multiple of 3")
            if len(r.seq) < start + 2:
                rep.errors.append(f"{tag}: too short for {TARGET_LABELS[key]} (needs >= {start + 2} nt)")
                continue
            internal = [i for i in range(0, len(r.seq) - 3, 3) if CODON_TABLE.get(r.seq[i:i + 3]) == "*"]
            if internal:
                rep.errors.append(f"{tag}: internal stop codon at codon {internal[0] // 3 + 1}")
            aa = translate(r.seq[start - 1:start + 2])
            if aa != EXPECTED_WT[key]:
                rep.warnings.append(
                    f"{tag}: {TARGET_LABELS[key]} reads {aa}, expected wild-type {EXPECTED_WT[key]}. "
                    "Check the reference is a susceptible wild-type CDS."
                )
    for gene in PMQR_GENES:
        for r in refs.get(gene, []):
            if len(r.seq) < 200:
                rep.warnings.append(f"{r.label}: unusually short {gene} reference ({len(r.seq)} nt)")
    return rep


# ------------------------------------------------------------------------------ engine
class SeptaTyperEngine:
    def __init__(
        self,
        reference_dir: str | os.PathLike = "references",
        blastn: str = "blastn",
        makeblastdb: str = "makeblastdb",
        threads: int = 2,
        timeout: int = 180,
        min_gene_identity: float = MIN_GENE_IDENTITY,
        max_bases: int = 30_000_000,
        max_contigs: int = 50_000,
    ):
        self.ref_dir = Path(reference_dir)
        self.blastn, self.makeblastdb = blastn, makeblastdb
        self.threads, self.timeout = threads, timeout
        self.min_gene_identity = min_gene_identity
        self.max_bases, self.max_contigs = max_bases, max_contigs

        for exe in (blastn, makeblastdb):
            if shutil.which(exe) is None:
                raise EngineError(f"'{exe}' not found on PATH. Install BLAST+ (e.g. apt install ncbi-blast+).")
        self.refs, rep = _load_refs(self.ref_dir)
        self.reference_report = validate_references(self.refs, rep)
        if not self.reference_report.ok:
            raise EngineError("Reference validation failed:\n  " + "\n  ".join(self.reference_report.errors))
        for w in self.reference_report.warnings:
            log.warning(w)

    # ---- helpers -------------------------------------------------------------------
    def _run(self, args: List[str]) -> str:
        try:
            p = subprocess.run(args, capture_output=True, text=True, timeout=self.timeout, check=False)
        except subprocess.TimeoutExpired as exc:
            raise EngineError(f"{Path(args[0]).name} timed out after {self.timeout}s") from exc
        if p.returncode != 0:
            raise EngineError(f"{Path(args[0]).name} failed: {p.stderr.strip()[:500]}")
        return p.stdout

    def _prepare_genome(self, text: str, workdir: Path) -> Tuple[Path, Dict[str, str], int, int]:
        records = read_fasta(text)
        if not records:
            raise InputError("No FASTA records found")
        if len(records) > self.max_contigs:
            raise InputError(f"Too many contigs ({len(records)}); limit is {self.max_contigs}")
        names: Dict[str, str] = {}
        lines: List[str] = []
        total = 0
        bad = 0
        for i, (header, seq) in enumerate(records, start=1):
            seq = re.sub(r"\s+", "", seq).upper()
            if not seq:
                continue
            bad += sum(1 for ch in seq if ch not in IUPAC_NT)
            total += len(seq)
            if total > self.max_bases:
                raise InputError(f"Assembly exceeds the {self.max_bases:,} base limit")
            cid = f"c{i}"
            names[cid] = header.split()[0] if header else cid
            lines.append(f">{cid}\n{seq}\n")
        if total == 0:
            raise InputError("All FASTA records are empty")
        if bad / total > 0.02:
            raise InputError("Sequence is not nucleotide FASTA (too many non-nucleotide characters)")
        path = workdir / "genome.fa"
        path.write_text("".join(lines))
        return path, names, len(names), total

    def _blast(self, query_records: List[Tuple[str, str]], db: Path, workdir: Path, tag: str) -> List[HSP]:
        q = workdir / f"{tag}.fa"
        q.write_text("".join(f">{i}\n{s}\n" for i, s in query_records))
        out = self._run([
            self.blastn, "-task", "blastn", "-query", str(q), "-db", str(db),
            "-evalue", "1e-10", "-perc_identity", "70", "-dust", "no", "-soft_masking", "false",
            "-max_target_seqs", "500", "-num_threads", str(self.threads),
            "-outfmt", f"6 {BLAST_FIELDS}",
        ])
        return parse_blast_tabular(out)

    # ---- chromosomal genes -----------------------------------------------------------
    def _candidate_refs(self, gene: str, species: Optional[str]) -> List[Tuple[int, Ref]]:
        allrefs = list(enumerate(self.refs[gene]))
        if species:
            same = [(i, r) for i, r in allrefs if r.species == species]
            if same:
                return same
        return allrefs

    def _analyse_genes(
        self, hsps: List[HSP], cand: Dict[str, List[Tuple[int, Ref]]], names: Dict[str, str], warnings: List[str]
    ) -> Tuple[Dict[str, dict], Dict[str, dict], Dict[str, Optional[Ref]]]:
        gene_info: Dict[str, dict] = {}
        codons: Dict[str, dict] = {}
        best_refs: Dict[str, Optional[Ref]] = {}
        gene_hsps: Dict[str, List[HSP]] = {}

        for gene in GENES:
            by_ref: Dict[int, List[HSP]] = defaultdict(list)
            for h in hsps:
                if h.gene == gene and h.pident >= self.min_gene_identity:
                    by_ref[h.ref_idx].append(h)
            if not by_ref:
                gene_info[gene] = {"found": False}
                best_refs[gene] = None
                gene_hsps[gene] = []
                warnings.append(f"{gene} was not found in the assembly (no hit >= {self.min_gene_identity:.0f}% identity)")
                continue
            best_idx = max(by_ref, key=lambda i: sum(h.bitscore for h in by_ref[i]))
            chosen = by_ref[best_idx]
            ref = self.refs[gene][best_idx]
            ident, cov = summarise_hsps(chosen)
            top = max(chosen, key=lambda h: h.bitscore)
            gene_info[gene] = {
                "found": True,
                "reference": ref.label,
                "reference_species": ref.species or ref.species_raw or None,
                "identity": ident,
                "coverage": cov,
                "contig": names.get(top.sseqid, top.sseqid),
                "start": min(top.sstart, top.send),
                "end": max(top.sstart, top.send),
                "strand": top.strand,
            }
            if ident < 95.0:
                warnings.append(f"{gene} is divergent from the closest reference ({ident:.1f}% identity)")
            best_refs[gene] = ref
            gene_hsps[gene] = chosen

        for key, (gene, codon_no, start) in TARGETS.items():
            label = TARGET_LABELS[key]
            call = {"position": label, "gene": gene, "ref_nt_start": start, "codon": None, "aa": None, "status": "gene_not_found"}
            candidates: List[Tuple[float, str, str]] = []
            status_seen = "not_covered"
            for h in gene_hsps.get(gene, []):
                codon, st = extract_codon(h, start)
                if codon:
                    candidates.append((h.bitscore, codon, h.sseqid))
                elif st != "not_covered":
                    status_seen = st
            if gene_hsps.get(gene):
                if candidates:
                    candidates.sort(reverse=True)
                    codon = candidates[0][1]
                    aa = translate(codon)
                    call.update(codon=codon, aa=aa, status="ok")
                    if aa == "*":
                        call.update(aa=None, status="stop_codon")
                        warnings.append(f"{label} is a stop codon ({codon})")
                    others = {c for _, c, _ in candidates[1:] if c != codon}
                    if others:
                        warnings.append(
                            f"{label}: multiple gene copies disagree ({codon} used; also {', '.join(sorted(others))})"
                        )
                else:
                    call["status"] = status_seen
                    warnings.append(f"{label} could not be called ({status_seen.replace('_', ' ')})")
            codons[key] = call
        return gene_info, codons, best_refs

    # ---- PMQR --------------------------------------------------------------------------
    def _call_pmqr(self, hsps: List[HSP], names: Dict[str, str], warnings: List[str]) -> Dict[str, dict]:
        out: Dict[str, dict] = {}
        for gene in PMQR_GENES:
            groups: Dict[Tuple[int, str, str], List[HSP]] = defaultdict(list)
            for h in hsps:
                if h.gene == gene:
                    groups[(h.ref_idx, h.sseqid, h.strand)].append(h)
            best = None
            for (ridx, sseqid, strand), hs in groups.items():
                ident, cov = summarise_hsps(hs)
                passes = ident >= PMQR_MIN_IDENTITY and cov >= PMQR_MIN_COVERAGE
                rank = (passes, ident * cov / 100.0)
                if best is None or rank > best[0]:
                    best = (rank, ident, cov, sseqid, ridx)
            if best is None:
                out[gene] = {"present": 0, "identity": None, "coverage": None, "contig": None, "closest_allele": None}
                continue
            (passes, _), ident, cov, sseqid, ridx = best
            out[gene] = {
                "present": int(passes), "identity": ident, "coverage": cov,
                "contig": names.get(sseqid, sseqid), "closest_allele": self.refs[gene][ridx].label,
            }
            if not passes and ident >= PMQR_MIN_IDENTITY and cov >= 50.0:
                warnings.append(f"{gene}: partial hit ({ident:.1f}% identity, {cov:.1f}% coverage) counted as absent")
        return out

    # ---- species -----------------------------------------------------------------------
    @staticmethod
    def _infer_species(best_refs: Mapping[str, Optional[Ref]], gene_info: Mapping[str, dict]) -> Tuple[Optional[str], str]:
        votes: Counter = Counter()
        n = 0
        for gene, ref in best_refs.items():
            if ref is not None and ref.species:
                votes[ref.species] += gene_info[gene].get("identity", 0.0)
                n += 1
        if not votes:
            return None, "0/3"
        sp = votes.most_common(1)[0][0]
        agree = sum(1 for r in best_refs.values() if r is not None and r.species == sp)
        return sp, f"{agree}/3"

    # ---- public API ----------------------------------------------------------------------
    def predict(self, fasta_text: str, species: Optional[str] = None) -> dict:
        user_species = canonical_species(species) if species else None
        if species and user_species is None:
            raise InputError(f"Unsupported species '{species}'. Supported: {', '.join(SPECIES)}")
        warnings: List[str] = []

        with tempfile.TemporaryDirectory(prefix="septatyper_") as tmp:
            work = Path(tmp)
            genome, names, n_contigs, total = self._prepare_genome(fasta_text, work)
            if total < 1_000_000:
                warnings.append(
                    f"Input is only {total:,} bases; a complete assembly is expected. Missing genes give an Indeterminate call."
                )
            db = work / "genomedb"
            self._run([self.makeblastdb, "-in", str(genome), "-dbtype", "nucl", "-out", str(db)])

            cand = {g: self._candidate_refs(g, None) for g in GENES}  # all references; species only drives the grammar
            chrom_q = [(f"{g}|{i}", r.seq) for g in GENES for i, r in cand[g]]
            chrom_hsps = self._blast(chrom_q, db, work, "chrom")
            gene_info, codons, best_refs = self._analyse_genes(chrom_hsps, cand, names, warnings)

            pmqr_q = [(f"{g}|{i}", r.seq) for g in PMQR_GENES for i, r in enumerate(self.refs[g])]
            pmqr_hsps = self._blast(pmqr_q, db, work, "pmqr")
            pmqr = self._call_pmqr(pmqr_hsps, names, warnings)

        # species
        inferred, agree = self._infer_species(best_refs, gene_info)
        if user_species:
            species_name, source = user_species, "user"
        else:
            species_name, source = inferred, ("inferred" if inferred else "unresolved")

        aa = {k: c["aa"] for k, c in codons.items()}
        recn_codon = codons["recN_346"]["codon"]
        recn = {
            "codon": recn_codon,
            "aa": aa["recN_346"],
            "leu_codon": recn_codon if aa["recN_346"] == "L" else None,
            "synonymous_flag": (recn_codon in ("CTA", "CTG")) if (species_name == SENTERICA and recn_codon) else None,
        }

        gr = run_grammar(aa, {g: pmqr[g]["present"] for g in PMQR_GENES}, species_name, recn_codon)
        if gr.prediction == "Indeterminate":
            warnings.append("Prediction is Indeterminate: one or more of GyrA-83, GyrA-87, ParC-80, ParC-84 could not be called")
        if aa["recN_346"] is None:
            warnings.append("RecN-346 could not be called, so RecN-dependent overrides were not evaluated")
        if source != "user" and species_name in (ECOLI, SSONNEI) and gr.tier1 and aa.get("recN_346") == "F" and \
                "-".join(str(aa[k]) for k in ("gyrA_83", "gyrA_87", "parC_80", "parC_84")) == "S-D-S-E":
            warnings.append(
                "RecN-346 F with S-D-S-E matches the S. sonnei override, but E. coli and S. sonnei cannot be told apart "
                "from gyrA/parC/recN alone. Resubmit with the species set to apply or exclude it."
            )

        return {
            "pipeline_version": PIPELINE_VERSION,
            "input": {"n_contigs": n_contigs, "total_bases": total},
            "species": {"name": species_name, "source": source, "genes_agreeing": agree, "inferred": inferred},
            "septatype": gr.septatype,
            "septatype_aa_level": gr.septatype_aa_level,
            "prediction": gr.prediction,
            "predicted_phenotype": gr.phenotype,
            "mean_hydrophobicity": None if gr.mean_hydrophobicity is None else round(gr.mean_hydrophobicity, 4),
            "hydrophobicity_scale": HYDROPHOBICITY_SCALE,
            "codons": codons,
            "recn346": recn,
            "pmqr": pmqr,
            "genes": gene_info,
            "tier1": None if gr.tier1 is None else {
                "baseline": gr.tier1.baseline, "votes": gr.tier1.votes,
                "decided_by": gr.tier1.decided_by, "positions": gr.tier1.detail,
            },
            "tier2": None if gr.tier2 is None else {
                "override_applied": gr.tier2.override_applied, "rules": gr.tier2.rules,
            },
            "warnings": warnings,
        }


# ---------------------------------------------------------------------------------- CLI
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="SeptaTyper prediction engine")
    ap.add_argument("genome", nargs="?", help="Genome assembly FASTA")
    ap.add_argument("--ref-dir", default=os.getenv("SEPTATYPER_REF_DIR", "references"))
    ap.add_argument("--species", default=None)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--check-refs", action="store_true", help="Validate the reference files and exit")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.check_refs:
        refs, rep = _load_refs(Path(args.ref_dir))
        rep = validate_references(refs, rep)
        for e in rep.errors:
            print("ERROR  ", e)
        for w in rep.warnings:
            print("WARNING", w)
        print("References OK" if rep.ok else "References FAILED")
        return 0 if rep.ok else 1
    if not args.genome:
        ap.error("genome FASTA required")
    engine = SeptaTyperEngine(args.ref_dir, threads=args.threads)
    result = engine.predict(Path(args.genome).read_text(), species=args.species)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
