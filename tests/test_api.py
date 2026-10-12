import gzip
import random

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from septatyper.database import get_session
from septatyper.models import Base, Septatype
from septatyper.routes import predict as predict_route
from tests.synthetic import variant, make_genome, plant


@pytest.fixture()
def client(refs, tmp_path, monkeypatch):
    monkeypatch.setenv("SEPTATYPER_REF_DIR", str(refs[0]))
    predict_route.get_prediction_engine.cache_clear()
    eng = create_engine(f"sqlite:///{tmp_path/'t.db'}", connect_args={"check_same_thread": False})
    Base.metadata.create_all(eng)
    Session = sessionmaker(eng)
    with Session() as s:
        s.add(Septatype(septatype_string="S-D-S-E-L-0-0", gyra83_aa="S", gyra87_aa="D", parc80_aa="S",
                        parc84_aa="E", recn346_aa="L", qnrb2=0, qnrs1=0, isolate_count=1129))
        s.commit()

    def override():
        with Session() as s:
            yield s

    from septatyper.main import app
    app.dependency_overrides[get_session] = override
    yield TestClient(app)
    app.dependency_overrides.clear()
    predict_route.get_prediction_engine.cache_clear()


@pytest.fixture()
def genome(refs):
    rng = random.Random(3)
    g = refs[1]
    seqs = {k: variant(g[k]["Klebsiella_pneumoniae"], 0.02, rng, k) for k in ("gyrA", "parC", "recN")}
    return make_genome(seqs, rng), seqs, rng


def post(client, text, **data):
    return client.post("/api/v1/predict", files={"file": ("g.fasta", text.encode())}, data=data)


def test_upload_wild_type_known_septatype(client, genome):
    r = post(client, genome[0], species="Klebsiella pneumoniae")
    assert r.status_code == 200, r.text
    j = r.json()
    assert j["septatype"] == "S-D-S-E-L-0-0" and j["prediction"] == "Csus"
    assert j["codons"]["gyrA_83"] == {"position": "GyrA-83", "gene": "gyrA", "ref_nt_start": 247,
                                      "codon": "TCG", "aa": "S", "status": "ok"}
    assert j["mean_hydrophobicity"] == pytest.approx(-0.5)
    assert j["database"] == {"in_database": True, "septatype_id": 1, "isolate_count": 1129, "is_novel": False}


def test_novel_septatype_is_flagged(client, genome):
    text, seqs, rng = genome
    mutated = make_genome({**seqs, "gyrA": plant(seqs["gyrA"], "gyrA_83", "ATC")}, rng)
    j = post(client, mutated, species="Klebsiella pneumoniae").json()
    assert j["septatype"] == "I-D-S-E-L-0-0" and j["database"]["is_novel"] is True


def test_gzip_and_pasted_text(client, genome):
    gz = gzip.compress(genome[0].encode())
    r = client.post("/api/v1/predict", files={"file": ("g.fa.gz", gz)})
    assert r.status_code == 200 and r.json()["prediction"] == "Csus"
    r2 = client.post("/api/v1/predict", data={"sequence": genome[0]})
    assert r2.status_code == 200


@pytest.mark.parametrize("kwargs,code", [
    ({}, 422),                                              # neither file nor sequence
    ({"data": {"sequence": "ACGT", "species": "Vibrio cholerae"}}, 422),
])
def test_validation_errors(client, kwargs, code):
    assert client.post("/api/v1/predict", **kwargs).status_code == code


def test_both_file_and_sequence(client):
    r = client.post("/api/v1/predict", files={"file": ("g.fa", b">a\nACGT\n")}, data={"sequence": ">b\nACGT"})
    assert r.status_code == 422


def test_protein_fasta_rejected(client):
    assert post(client, ">p\n" + "MKVLAAGIVGLLLAQPSTAWQERTY" * 400).status_code == 400


def test_oversize_rejected(client, monkeypatch):
    monkeypatch.setattr(predict_route, "MAX_UPLOAD_BYTES", 1000)
    assert post(client, ">a\n" + "ACGT" * 1000).status_code == 413


def test_engine_unavailable_returns_503(client, monkeypatch, tmp_path):
    monkeypatch.setenv("SEPTATYPER_REF_DIR", str(tmp_path / "nope"))
    predict_route.get_prediction_engine.cache_clear()
    r = post(client, ">a\nACGT\n")
    assert r.status_code == 503 and "Missing reference file" in r.json()["detail"]
