import pytest

from septatyper import prediction_engine as pe
from septatyper.prediction_engine import ECOLI, KPNEU, SENTERICA, SSONNEI, run_grammar


def aa(t, recn="L"):
    g83, g87, p80, p84 = t.split("-")
    return {"gyrA_83": g83, "gyrA_87": g87, "parC_80": p80, "parC_84": p84, "recN_346": recn}


def go(t, qb=0, qs=0, sp=ECOLI, recn="L", codon="CTG"):
    return run_grammar(aa(t, recn), {"qnrB2": qb, "qnrS1": qs}, sp, codon)


@pytest.mark.parametrize("tet,expected", [
    ("S-D-S-E", "Csus"),    # fully polar/negative
    ("I-D-I-E", "Cres"),
    ("I-N-I-E", "Cres"),
    ("S-D-I-E", "Cres"),    # ParC-84 E paired with hydrophobic ParC-80 (I-E) is resistant
    ("F-N-S-K", "Cres"),
    ("L-D-S-E", "Csus"),    # one vote, <H> below the scale mean
    ("Y-D-S-E", "Csus"),    # Y counts as polar
    ("S-N-S-E", "Csus"),
    ("L-N-I-V", "Cres"),
])
def test_tier1(tet, expected):
    assert go(tet).prediction == expected


def test_septatype_string_and_hydrophobicity():
    r = go("S-D-S-E", qs=1, sp=KPNEU)
    assert r.septatype == "S-D-S-E-L-0-1"
    assert r.mean_hydrophobicity == pytest.approx((-0.18 - 0.90 - 0.18 - 0.74) / 4)  # = -0.5


def test_ParC84_E_is_context_dependent():
    s_e = go("S-D-S-E").tier1
    i_e = go("S-D-I-E").tier1
    assert s_e.votes == 0 and i_e.votes == 2


@pytest.mark.parametrize("qb,qs", [(1, 0), (0, 1), (1, 1)])
def test_qnr_override(qb, qs):
    r = go("S-D-S-E", qb=qb, qs=qs)
    assert r.prediction == "Cres" and r.tier2.rules == ["qnr_present"]


def test_sonnei_recN_override():
    assert go("S-D-S-E", sp=SSONNEI, recn="F").tier2.rules == ["sonnei_recN_L346F"]
    assert go("S-D-S-E", sp=SSONNEI, recn="L").prediction == "Csus"
    assert go("S-D-S-E", sp=ECOLI, recn="F").prediction == "Csus"


def test_enterica_ctg_override():
    assert go("F-D-S-E", sp=SENTERICA, codon="CTG").tier2.rules == ["enterica_recN_ctg"]
    assert go("F-D-S-E", sp=SENTERICA, codon="CTA").prediction == "Csus"
    assert go("F-D-S-E", sp=ECOLI, codon="CTG").prediction == "Csus"


def test_override_never_downgrades():
    r = go("I-N-I-E", sp=SENTERICA)
    assert r.prediction == "Cres" and not r.tier2.override_applied


def test_hydrophobicity_tiebreak(monkeypatch):
    monkeypatch.setattr(pe, "H_RESCUE_THRESHOLD", -0.3)
    r = go("L-D-S-E")  # <H> = (1.06 - 0.90 - 0.18 - 0.74) / 4 = -0.19
    assert r.prediction == "Cres" and r.tier1.decided_by == "hydrophobicity"


def test_missing_residue_is_indeterminate():
    a = aa("S-D-S-E"); a["parC_80"] = None
    r = run_grammar(a, {"qnrB2": 0, "qnrS1": 0}, ECOLI, "CTG")
    assert r.prediction == "Indeterminate" and r.septatype is None


@pytest.mark.parametrize("tet", ["L-D-S-E", "S-Y-S-E"])
def test_sonnei_single_vote_override(tet):
    r = go(tet, sp=SSONNEI)
    assert r.prediction == "Cres" and r.tier2.rules == ["sonnei_single_vote_override"]
    assert go(tet, sp=ECOLI).prediction == "Csus"          # species-specific
    assert go(tet, sp=KPNEU).prediction == "Csus"


def test_sonnei_override_scope():
    assert go("S-N-S-E", sp=SSONNEI).prediction == "Csus"   # other single-vote tetratypes unaffected
    assert go("S-D-S-E", sp=SSONNEI).prediction == "Csus"
    r = go("L-D-S-E", qs=1, sp=SSONNEI)                     # rules accumulate; outcome stays Resistant
    assert r.prediction == "Cres" and r.tier2.rules == ["qnr_present", "sonnei_single_vote_override"]


def test_codon_notation_only_for_enterica_and_klebsiella_FDSE():
    assert go("F-D-S-E", sp=SENTERICA, codon="CTG").septatype == "F-D-S-E-Lctg-0-0"
    assert go("F-D-S-E", sp=SENTERICA, codon="CTA").septatype == "F-D-S-E-Lcta-0-0"
    assert go("F-D-S-E", sp=SENTERICA, codon="CTG").septatype_aa_level == "F-D-S-E-L-0-0"
    assert go("F-D-S-E", qs=1, sp=SENTERICA, codon="CTG").septatype == "F-D-S-E-Lctg-0-1"
    assert go("F-D-S-E", sp=SENTERICA, codon="CTT").septatype == "F-D-S-E-L-0-0"      # other Leu codons stay plain
    assert go("F-D-S-E", sp=KPNEU, codon="CTA").septatype == "F-D-S-E-Lcta-0-0"        # K. pneumoniae uses it too
    assert go("F-D-S-E", qs=1, sp=KPNEU, codon="CTG").septatype == "F-D-S-E-Lctg-0-1"
    assert go("F-D-S-E", sp=ECOLI, codon="CTA").septatype == "F-D-S-E-L-0-0"          # other species stay plain
    assert go("F-D-S-E", sp=KPNEU, codon="CTG").prediction == "Csus"                  # notation only, no new override
    assert go("S-D-S-E", sp=SENTERICA, codon="CTG").septatype == "S-D-S-E-L-0-0"      # other tetratypes stay plain
    assert go("F-D-S-E", sp=SENTERICA, recn="F", codon="TTC").septatype == "F-D-S-E-F-0-0"


def test_eisenberg_values_match_manuscript():
    assert (pe.EISENBERG["S"], pe.EISENBERG["I"], pe.EISENBERG["L"], pe.EISENBERG["F"]) == (-0.18, 1.38, 1.06, 1.19)
    assert abs(sum(pe.EISENBERG.values())) < 1e-9          # zero-mean scale, so the 0.0 tiebreak is the scale average
