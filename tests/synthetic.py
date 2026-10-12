"""Synthetic references and genomes for testing the engine without real sequence data.

These sequences are random. They exist only to exercise coordinates, strand handling,
indels, contig breaks and thresholds. Real references must be supplied for production.
"""
import random
from pathlib import Path
from typing import Dict, Optional

from septatyper.prediction_engine import CODON_TABLE, TARGETS

SENSE = [c for c, a in CODON_TABLE.items() if a != "*"]
LENGTHS = {"gyrA": 2628, "parC": 2253, "recN": 1668}
WT = {"gyrA_83": "TCG", "gyrA_87": "GAC", "parC_80": "AGC", "parC_84": "GAA", "recN_346": "CTG"}
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def random_dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_gene(gene: str, rng: random.Random) -> str:
    codons = ["ATG"] + [rng.choice(SENSE) for _ in range(LENGTHS[gene] // 3 - 1)]
    seq = "".join(codons)
    for key, (g, _, start) in TARGETS.items():
        if g == gene:
            seq = plant(seq, key, WT[key])
    return seq


def plant(seq: str, key: str, codon: str) -> str:
    start = TARGETS[key][2]
    return seq[: start - 1] + codon + seq[start + 2:]


def variant(seq: str, rate: float, rng: random.Random, gene: str) -> str:
    """Point-mutate outside the +-30 nt windows around target codons; never create a stop."""
    protected = [(s - 30, s + 32) for g, _, s in TARGETS.values() if g == gene]
    out = list(seq)
    for i in range(3, len(seq)):
        if any(a <= i + 1 <= b for a, b in protected) or rng.random() > rate:
            continue
        old = out[i]
        out[i] = rng.choice([b for b in "ACGT" if b != old])
        c0 = (i // 3) * 3
        if CODON_TABLE["".join(out[c0:c0 + 3])] == "*":
            out[i] = old
    return "".join(out)


def make_reference_dir(path: Path, seed: int = 7) -> Dict[str, Dict[str, str]]:
    rng = random.Random(seed)
    (path / "pmqr").mkdir(parents=True, exist_ok=True)
    species_div = {
        "Escherichia_coli": 0.0, "Shigella_sonnei": 0.004,
        "Salmonella_enterica": 0.05, "Klebsiella_pneumoniae": 0.08,
    }
    genes: Dict[str, Dict[str, str]] = {}
    for gene in LENGTHS:
        base = make_gene(gene, rng)
        genes[gene] = {sp: (base if d == 0 else variant(base, d, rng, gene)) for sp, d in species_div.items()}
        with open(path / f"{gene}.fasta", "w") as fh:
            for sp, seq in genes[gene].items():
                fh.write(f">{gene}|{sp}|SYN_{sp}\n{seq}\n")
    pmqr = {"qnrB2": random_dna(645, rng), "qnrS1": random_dna(657, rng)}
    for g, s in pmqr.items():
        (path / "pmqr" / f"{g}.fasta").write_text(f">{g}|ref|SYN\n{s}\n")
    genes["_pmqr"] = pmqr
    return genes


def diverge(seq: str, identity: float, rng: random.Random) -> str:
    out = list(seq)
    for i in rng.sample(range(len(seq)), int(len(seq) * (1 - identity))):
        out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def make_genome(
    gene_seqs: Dict[str, str],
    rng: random.Random,
    minus: Optional[set] = None,
    split: Optional[Dict[str, int]] = None,
    extra: Optional[Dict[str, str]] = None,
    flank: int = 20000,
) -> str:
    """Build a FASTA assembly. `minus` genes are reverse-complemented; `split` maps a gene
    to the 1-based nt where it is cut across two contigs; `extra` adds more sequences
    (e.g. qnr genes) on their own contigs."""
    minus = minus or set()
    split = split or {}
    contigs = []
    for gene, seq in gene_seqs.items():
        s = revcomp(seq) if gene in minus else seq
        if gene in split:
            cut = split[gene]
            cut = len(seq) - cut if gene in minus else cut
            contigs.append(random_dna(flank, rng) + s[:cut])
            contigs.append(s[cut:] + random_dna(flank, rng))
        else:
            contigs.append(random_dna(flank, rng) + s + random_dna(flank, rng))
    for name, seq in (extra or {}).items():
        contigs.append(random_dna(flank, rng) + seq + random_dna(flank, rng))
    contigs.append(random_dna(50000, rng))
    return "".join(f">contig_{i + 1} test\n{c}\n" for i, c in enumerate(contigs))
