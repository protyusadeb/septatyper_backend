import random

import pytest
import requests
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker

from septatyper.auth import hash_password
from septatyper.bvbrc import BVBRCClient
from septatyper.database import get_session
from septatyper.engine_runtime import get_prediction_engine
from septatyper.models import (Base, Isolate, IsolateGenotype, PMQRCall, Septatype, Submission, User, UserRole)
from septatyper.ratelimit import RateLimiter
from septatyper.routes import submissions as sub_route
from tests.synthetic import make_genome, plant, variant

META = dict(contact_email="Lab@Example.org", species="Klebsiella pneumoniae", habitat="Urine",
            country="Bangladesh", collection_year="2019", lab_phenotype="Susceptible")


@pytest.fixture()
def env(refs, tmp_path, monkeypatch):
    monkeypatch.setenv("SEPTATYPER_REF_DIR", str(refs[0]))
    monkeypatch.setenv("SEPTATYPER_JWT_SECRET", "x" * 40)
    get_prediction_engine.cache_clear()
    monkeypatch.setattr(sub_route, "SUBMISSION_DIR", tmp_path / "subs")
    monkeypatch.setattr(sub_route, "submit_limiter", RateLimiter(1000, 3600))
    eng = create_engine(f"sqlite:///{tmp_path/'t.db'}", connect_args={"check_same_thread": False})
    Base.metadata.create_all(eng)
    Session = sessionmaker(eng)
    with Session() as s:
        s.add(Septatype(septatype_string="S-D-S-E-L-0-0", gyra83_aa="S", gyra87_aa="D", parc80_aa="S",
                        parc84_aa="E", recn346_aa="L", qnrb2=0, qnrs1=0, isolate_count=0))
        s.add(User(email="cur@x.org", hashed_password=hash_password("curator-pass"), role=UserRole.CURATOR))
        s.add(User(email="sub@x.org", hashed_password=hash_password("submit-pass"), role=UserRole.SUBMITTER))
        s.commit()

    def override():
        with Session() as s:
            yield s

    from septatyper.main import app
    app.dependency_overrides[get_session] = override
    c = TestClient(app)
    yield c, Session
    app.dependency_overrides.clear()
    get_prediction_engine.cache_clear()


@pytest.fixture()
def bvbrc(monkeypatch):
    calls = []

    def fake(self, ids):
        calls.append(ids)
        return {i: {"genome_name": f"Klebsiella pneumoniae strain {i}",
                    "amr": [{"evidence": "Computational Method", "resistant_phenotype": "Resistant",
                             "computational_method": "AdaBoost", "computational_method_version": "1"}]} for i in ids}
    monkeypatch.setattr(BVBRCClient, "fetch_batch", fake)
    return calls


def genome(refs, muts=None, src="Klebsiella_pneumoniae", seed=5):
    rng = random.Random(seed)
    seqs = {}
    for g in ("gyrA", "parC", "recN"):
        s = variant(refs[1][g][src], 0.02, rng, g)
        for k, codon in (muts or {}).items():
            from septatyper.prediction_engine import TARGETS
            if TARGETS[k][0] == g:
                s = plant(s, k, codon)
        seqs[g] = s
    return make_genome(seqs, rng)


def submit(c, fasta, **over):
    data = {**META, **over}
    data = {k: v for k, v in data.items() if v is not None}
    return c.post("/api/v1/submissions", files={"file": ("my genome.fa", fasta.encode())}, data=data)


def login(c, email="cur@x.org", pw="curator-pass"):
    r = c.post("/api/v1/auth/login", json={"email": email, "password": pw})
    return {"Authorization": f"Bearer {r.json()['access_token']}"} if r.status_code == 200 else None


def test_full_flow_clean_submission_then_approve(env, refs, bvbrc):
    c, Session = env
    r = submit(c, genome(refs), bv_brc_id="562.9999", isolate_name="BD-17", mic_value="0.06")
    assert r.status_code == 202, r.text
    sid, token = r.json()["submission_id"], r.json()["status_token"]
    assert bvbrc == [["562.9999"]]

    st = c.get(f"/api/v1/submissions/{sid}/status", params={"token": token}).json()
    assert st["analysis_state"] == "done" and st["status"] == "pending" and st["flags"] == []
    assert st["bv_brc_name_fetched"] == "Klebsiella pneumoniae strain 562.9999"
    assert st["septatype"] == "S-D-S-E-L-0-0" and st["septatyper_prediction"] == "Csus"

    # staged only: the live tables are untouched
    with Session() as s:
        assert s.scalar(select(Isolate)) is None

    h = login(c)
    q = c.get("/api/v1/submissions/pending", headers=h).json()
    assert q["total"] == 1 and q["counts"] == {"needs_review": 0, "pending": 1, "failed": 0}
    d = c.get(f"/api/v1/submissions/{sid}", headers=h).json()
    assert d["contact_email"] == "lab@example.org" and d["prediction"]["codons"]["gyrA_83"]["codon"] == "TCG"

    r = c.post(f"/api/v1/submissions/{sid}/approve", headers=h, json={"notes": "checked"})
    assert r.status_code == 200 and r.json()["isolate_id"] == "562.9999" and r.json()["created_septatype"] is False
    with Session() as s:
        iso = s.get(Isolate, "562.9999")
        assert iso.source.value == "community" and iso.official_genome_name.endswith("562.9999")
        assert iso.lab_phenotype.value == "Susceptible" and iso.septatyper_prediction.value == "Csus"
        assert iso.bv_brc_computational_prediction.value == "Resistant"
        assert iso.septatyper_wins is True and iso.bvbrc_correct is False   # SeptaTyper right, BV-BRC wrong
        assert iso.country_iso3 == "BGD" and iso.collection_year_start == 2019 and float(iso.mic_value) == 0.06
        g = s.get(IsolateGenotype, "562.9999")
        assert (g.gyra83_codon, g.gyra83_aa, g.recn346_codon) == ("TCG", "S", "CTG")
        assert len(s.scalars(select(PMQRCall).where(PMQRCall.bv_brc_genome_id == "562.9999")).all()) == 2
        assert s.scalar(select(Septatype.isolate_count)) == 1
        assert s.get(Submission, sid).status.value == "approved"
    assert c.post(f"/api/v1/submissions/{sid}/approve", headers=h, json={}).status_code == 409
    assert c.get("/api/v1/isolates/562.9999").json()["septatype"] == "S-D-S-E-L-0-0"


def test_novel_and_mismatch_flags_require_confirmation(env, refs, bvbrc):
    c, Session = env
    r = submit(c, genome(refs, {"gyrA_83": "ATC"}), lab_phenotype="Resistant")
    sid, token = r.json()["submission_id"], r.json()["status_token"]
    st = c.get(f"/api/v1/submissions/{sid}/status", params={"token": token}).json()
    codes = {f["code"] for f in st["flags"]}
    assert st["status"] == "needs_review" and st["septatype"] == "I-D-S-E-L-0-0"
    assert {"novel_septatype", "phenotype_mismatch"} <= codes
    h = login(c)
    q = c.get("/api/v1/submissions/pending", params={"novel": "true", "mismatch": "true"}, headers=h).json()
    assert q["total"] == 1 and q["items"][0]["is_novel_septatype"] is True
    assert q["items"][0]["prediction_matches_phenotype"] is False

    no = c.post(f"/api/v1/submissions/{sid}/approve", headers=h, json={})
    assert no.status_code == 409 and "confirm" in no.json()["detail"]
    ok = c.post(f"/api/v1/submissions/{sid}/approve", headers=h, json={"confirm_flags": True})
    assert ok.status_code == 200 and ok.json()["created_septatype"] is True
    iid = ok.json()["isolate_id"]
    assert iid == f"SEPTA-{sid:06d}"                       # no BV-BRC ID supplied
    with Session() as s:
        sept = s.scalar(select(Septatype).where(Septatype.septatype_string == "I-D-S-E-L-0-0"))
        assert sept.is_novel is True and sept.isolate_count == 1 and sept.mean_hydrophobicity is not None
    assert c.get(f"/api/v1/isolates/{iid}").json()["bv_brc_url"] is None


def test_bvbrc_not_found_and_service_error(env, refs, monkeypatch):
    c, _ = env
    monkeypatch.setattr(BVBRCClient, "fetch_batch", lambda self, ids: {i: {"genome_name": None, "amr": []} for i in ids})
    r = submit(c, genome(refs), bv_brc_id="1.2")
    st = c.get(f"/api/v1/submissions/{r.json()['submission_id']}/status", params={"token": r.json()["status_token"]}).json()
    assert "bvbrc_id_not_found" in {f["code"] for f in st["flags"]} and st["bv_brc_name_fetched"] is None

    def boom(self, ids):
        raise requests.ConnectionError("down")
    monkeypatch.setattr(BVBRCClient, "fetch_batch", boom)
    r2 = submit(c, genome(refs, seed=6), bv_brc_id="3.4")
    st2 = c.get(f"/api/v1/submissions/{r2.json()['submission_id']}/status", params={"token": r2.json()["status_token"]}).json()
    assert st2["analysis_state"] == "done" and "bvbrc_lookup_failed" in {f["code"] for f in st2["flags"]}


def test_species_mismatch_duplicate_and_mic_conflict(env, refs, bvbrc):
    c, _ = env
    r = submit(c, genome(refs, src="Salmonella_enterica"), mic_value="2", mic_unit="mg/L")  # declared Klebsiella
    codes = {f["code"] for f in c.get(f"/api/v1/submissions/{r.json()['submission_id']}/status",
                                      params={"token": r.json()["status_token"]}).json()["flags"]}
    assert {"species_mismatch", "mic_phenotype_conflict"} <= codes
    r2 = submit(c, genome(refs, src="Salmonella_enterica"))
    codes2 = {f["code"] for f in c.get(f"/api/v1/submissions/{r2.json()['submission_id']}/status",
                                       params={"token": r2.json()["status_token"]}).json()["flags"]}
    assert "duplicate_genome" in codes2


def test_indeterminate_cannot_be_approved(env, bvbrc):
    c, _ = env
    rng = random.Random(1)
    junk = make_genome({}, rng)                                  # random contig, no target genes
    r = submit(c, junk)
    sid = r.json()["submission_id"]
    d = c.get(f"/api/v1/submissions/{sid}", headers=login(c)).json()
    assert "indeterminate_prediction" in {f["code"] for f in d["flags"]}
    assert c.post(f"/api/v1/submissions/{sid}/approve", headers=login(c), json={"confirm_flags": True}).status_code == 409


def test_reject_flow(env, refs, bvbrc):
    c, Session = env
    sid = submit(c, genome(refs)).json()["submission_id"]
    h = login(c)
    assert c.post(f"/api/v1/submissions/{sid}/reject", headers=h, json={"notes": "x"}).status_code == 422
    assert c.post(f"/api/v1/submissions/{sid}/reject", headers=h, json={"notes": "Low quality assembly"}).status_code == 204
    assert c.get("/api/v1/submissions/pending", headers=h).json()["total"] == 0
    assert c.post(f"/api/v1/submissions/{sid}/approve", headers=h, json={}).status_code == 409
    with Session() as s:
        assert s.scalar(select(Isolate)) is None


def test_access_control(env, refs, bvbrc):
    c, _ = env
    sid = submit(c, genome(refs)).json()["submission_id"]
    for method, url in [("get", "/api/v1/submissions/pending"), ("get", f"/api/v1/submissions/{sid}"),
                        ("post", f"/api/v1/submissions/{sid}/approve")]:
        assert getattr(c, method)(url).status_code == 401
    sub_h = login(c, "sub@x.org", "submit-pass")
    assert c.get("/api/v1/submissions/pending", headers=sub_h).status_code == 403
    assert c.post(f"/api/v1/submissions/{sid}/approve", headers=sub_h, json={}).status_code == 403
    assert c.post("/api/v1/auth/login", json={"email": "cur@x.org", "password": "wrong"}).status_code == 401
    assert c.post("/api/v1/auth/login", json={"email": "nobody@x.org", "password": "x"}).status_code == 401
    assert c.get("/api/v1/auth/me", headers={"Authorization": "Bearer garbage"}).status_code == 401
    assert c.get("/api/v1/auth/me", headers=login(c)).json() == {"email": "cur@x.org", "role": "curator"}
    assert c.get(f"/api/v1/submissions/{sid}/status", params={"token": "wrong-token-value"}).status_code == 404


@pytest.mark.parametrize("over,code", [
    ({"contact_email": "nope"}, 422), ({"collection_year": "1850"}, 422), ({"collection_year": "2999"}, 422),
    ({"country": "Atlantis"}, 422), ({"mic_value": "abc"}, 422), ({"mic_value": "-1"}, 422),
    ({"mic_value": "1", "mic_unit": "g/L"}, 422), ({"bv_brc_id": "abc"}, 422),
    ({"species": "Vibrio cholerae"}, 422), ({"lab_phenotype": "Maybe"}, 422),
])
def test_field_validation(env, refs, bvbrc, over, code):
    c, _ = env
    assert submit(c, genome(refs), **over).status_code == code


def test_file_validation_and_duplicates(env, refs, bvbrc, monkeypatch):
    c, Session = env
    assert submit(c, ">p\n" + "MKVLAAGIVGLLLAQPSTAWQERTY" * 40).status_code == 400
    assert c.post("/api/v1/submissions", data=META).status_code == 422                 # no file
    monkeypatch.setattr(sub_route, "MAX_UPLOAD_BYTES", 500)
    assert submit(c, ">a\n" + "ACGT" * 500).status_code == 413
    monkeypatch.setattr(sub_route, "MAX_UPLOAD_BYTES", 60 * 1024 * 1024)
    ok = submit(c, genome(refs), bv_brc_id="562.5")
    assert ok.status_code == 202
    assert submit(c, genome(refs, seed=9), bv_brc_id="562.5").status_code == 409        # already submitted
    from septatyper.models import Phenotype, Species
    with Session() as s:                                                                  # already live
        sp = Species(scientific_name="Escherichia coli")
        s.add(sp); s.flush()
        s.add(Isolate(bv_brc_genome_id="777.1", species_id=sp.species_id, lab_phenotype=Phenotype.RESISTANT))
        s.commit()
    assert submit(c, genome(refs, seed=10), bv_brc_id="777.1").status_code == 409


def test_rate_limit(env, refs, bvbrc, monkeypatch):
    c, _ = env
    monkeypatch.setattr(sub_route, "submit_limiter", RateLimiter(1, 3600))
    g = genome(refs)
    assert submit(c, g).status_code == 202
    assert submit(c, g).status_code == 429
