import pytest

from septatyper.prediction_engine import (
    EngineError, InputError, SeptaTyperEngine, TARGETS, _load_refs, validate_references,
)
from tests.synthetic import WT, diverge, make_genome, plant, random_dna, variant


def build(refs, rng, muts=None, species_src="Salmonella_enterica", rate=0.02, **kw):
    """Genome whose genes are diverged copies of the reference species, with codons planted."""
    _, genes = refs
    muts = muts or {}
    seqs = {}
    for gene in ("gyrA", "parC", "recN"):
        s = variant(genes[gene][species_src], rate, rng, gene)
        for key, codon in muts.items():
            if TARGETS[key][0] == gene:
                s = plant(s, key, codon)
        seqs[gene] = s
    return make_genome(seqs, rng, **kw)


def test_wild_type_forward(engine, refs, rng):
    r = engine.predict(build(refs, rng))
    assert r["septatype"] == "S-D-S-E-L-0-0" and r["prediction"] == "Csus"
    assert {k: c["codon"] for k, c in r["codons"].items()} == WT
    assert r["species"]["name"] == "Salmonella enterica subspecies I" and r["species"]["source"] == "inferred"


def test_known_mutations_and_reverse_strand(engine, refs, rng):
    muts = {"gyrA_83": "ATC", "gyrA_87": "AAC", "parC_80": "ATC"}
    fwd = engine.predict(build(refs, rng, muts))
    rev = engine.predict(build(refs, rng, muts, minus={"gyrA", "parC", "recN"}))
    for r in (fwd, rev):
        assert r["septatype"] == "I-N-I-E-L-0-0" and r["prediction"] == "Cres"
        assert r["codons"]["gyrA_83"]["codon"] == "ATC" and r["codons"]["parC_80"]["codon"] == "ATC"
    assert rev["genes"]["gyrA"]["strand"] == "-" and fwd["genes"]["gyrA"]["strand"] == "+"


def test_gene_split_across_contigs(engine, refs, rng):
    # cut gyrA between codon 83 (247-249) and codon 87 (259-261), parC between 80 and 84
    r = engine.predict(build(refs, rng, {"gyrA_87": "AAC"}, split={"gyrA": 255, "parC": 245}))
    assert r["codons"]["gyrA_83"]["codon"] == "TCG" and r["codons"]["gyrA_87"]["codon"] == "AAC"
    assert r["codons"]["parC_80"]["status"] == "ok" and r["codons"]["parC_84"]["status"] == "ok"
    r2 = engine.predict(build(refs, rng, split={"gyrA": 260}))  # cut inside codon 87
    assert r2["codons"]["gyrA_87"]["codon"] is None
    assert r2["prediction"] == "Indeterminate" and r2["septatype"] is None


def test_deletion_inside_codon_is_not_called(engine, refs, rng):
    _, genes = refs
    s = variant(genes["gyrA"]["Salmonella_enterica"], 0.02, rng, "gyrA")
    s = s[:247] + s[248:]  # delete the 2nd nt of codon 83
    others = {g: variant(genes[g]["Salmonella_enterica"], 0.02, rng, g) for g in ("parC", "recN")}
    r = engine.predict(make_genome({"gyrA": s, **others}, rng))
    assert r["codons"]["gyrA_83"]["codon"] is None and r["prediction"] == "Indeterminate"


def test_missing_gene(engine, refs, rng):
    _, genes = refs
    seqs = {g: variant(genes[g]["Salmonella_enterica"], 0.02, rng, g) for g in ("gyrA", "recN")}
    r = engine.predict(make_genome(seqs, rng))
    assert r["genes"]["parC"]["found"] is False and r["prediction"] == "Indeterminate"
    assert any("parC" in w for w in r["warnings"])


def test_pmqr_thresholds(engine, refs, rng):
    pm = refs[1]["_pmqr"]
    ok = diverge(pm["qnrB2"], 0.85, rng)                  # 85% identity, full length  -> present
    low = diverge(pm["qnrS1"], 0.72, rng)                  # 72% identity               -> absent
    r = engine.predict(build(refs, rng, extra={"qnrB2": ok, "qnrS1": low}))
    assert r["pmqr"]["qnrB2"]["present"] == 1 and r["pmqr"]["qnrS1"]["present"] == 0
    assert r["septatype"].endswith("-1-0")
    assert r["prediction"] == "Cres" and r["tier2"]["rules"] == ["qnr_present"]

    trunc = pm["qnrS1"][: int(len(pm["qnrS1"]) * 0.8)]     # 80% coverage at 100% identity -> absent
    r2 = engine.predict(build(refs, rng, extra={"qnrS1": trunc}))
    assert r2["pmqr"]["qnrS1"]["present"] == 0 and 79 <= r2["pmqr"]["qnrS1"]["coverage"] <= 81
    assert any("partial hit" in w for w in r2["warnings"])

    exact = engine.predict(build(refs, rng, extra={"qnrS1": pm["qnrS1"]}))
    assert exact["pmqr"]["qnrS1"]["present"] == 1


def test_sonnei_override_needs_species(engine, refs, rng):
    muts = {"recN_346": "TTC"}  # L -> F
    g = build(refs, rng, muts, species_src="Escherichia_coli")
    r_user = engine.predict(g, species="Shigella sonnei")
    assert r_user["prediction"] == "Cres" and r_user["tier2"]["rules"] == ["sonnei_recN_L346F"]
    r_ec = engine.predict(g, species="Escherichia coli")
    assert r_ec["prediction"] == "Csus"
    r_inf = engine.predict(g)  # inferred E. coli/S. sonnei: must warn, not guess
    assert any("S. sonnei override" in w for w in r_inf["warnings"])


def test_enterica_ctg_flag(engine, refs, rng):
    base = {"gyrA_83": "TTC"}  # F-D-S-E
    ctg = engine.predict(build(refs, rng, {**base, "recN_346": "CTG"}), species="Salmonella enterica subspecies I")
    cta = engine.predict(build(refs, rng, {**base, "recN_346": "CTA"}), species="Salmonella enterica subspecies I")
    assert ctg["recn346"]["synonymous_flag"] is True and ctg["prediction"] == "Cres"
    assert cta["recn346"]["synonymous_flag"] is True and cta["recn346"]["codon"] == "CTA" and cta["prediction"] == "Csus"
    assert cta["recn346"]["aa"] == "L" == ctg["recn346"]["aa"]
    assert ctg["septatype"] == "F-D-S-E-Lctg-0-0" and cta["septatype"] == "F-D-S-E-Lcta-0-0"
    assert cta["septatype_aa_level"] == "F-D-S-E-L-0-0"


def test_duplicate_copy_disagreement_warns(engine, refs, rng):
    _, genes = refs
    wt = variant(genes["gyrA"]["Salmonella_enterica"], 0.02, rng, "gyrA")
    mutant = plant(wt, "gyrA_83", "ATC")
    seqs = {g: variant(genes[g]["Salmonella_enterica"], 0.02, rng, g) for g in ("parC", "recN")}
    fa = make_genome({"gyrA": wt, **seqs}, rng, extra={"gyrA2": mutant})
    r = engine.predict(fa)
    assert any("multiple gene copies disagree" in w for w in r["warnings"])


def test_bad_inputs(engine):
    with pytest.raises(InputError):
        engine.predict("MKVLAAGIVGLLLAQPSTA\nMKVLAAGIVG")                      # no header
    with pytest.raises(InputError):
        engine.predict(">prot\n" + "MKVLAAGIVGLLLAQPSTAWQERTY" * 400)          # protein
    with pytest.raises(InputError):
        engine.predict(">x\nACGT" * 10, species="Vibrio cholerae")


def test_reference_validation(tmp_path, refs):
    import shutil
    d = tmp_path / "refs"
    shutil.copytree(refs[0], d)
    text = (d / "gyrA.fasta").read_text().replace("|SYN_Escherichia_coli\nATG", "|SYN_Escherichia_coli\nCCC", 1)
    (d / "gyrA.fasta").write_text(text)
    with pytest.raises(EngineError, match="does not start with ATG"):
        SeptaTyperEngine(d)
    rep = validate_references(*_load_refs(d)[:1])
    assert not rep.ok


def test_missing_blast_binary(refs):
    with pytest.raises(EngineError, match="not found on PATH"):
        SeptaTyperEngine(refs[0], blastn="definitely-not-blastn")
