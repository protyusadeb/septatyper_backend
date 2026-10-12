"""Fill isolates.septatyper_prediction for every isolate from its stored amino acids.

    python -m scripts.backfill_predictions [--dry-run]

Run after ingestion (and after any grammar change). The generated columns septatyper_correct and
septatyper_wins update automatically. Table S1 has no codons, so the S. enterica CTG rule cannot fire
for ingested isolates until codons are filled by the Module A pipeline.
"""
import argparse
from collections import Counter

from sqlalchemy import select

from septatyper.database import SessionLocal, get_engine
from septatyper.models import Isolate, IsolateGenotype, PMQRCall, SeptaPrediction, Species
from septatyper.prediction_engine import run_grammar


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    get_engine()
    stats: Counter = Counter()
    with SessionLocal() as db:
        pm = {}
        for c in db.scalars(select(PMQRCall)):
            pm.setdefault(c.bv_brc_genome_id, {})[c.gene] = int(c.present)
        rows = db.execute(
            select(Isolate, IsolateGenotype, Species.scientific_name)
            .join(IsolateGenotype, IsolateGenotype.bv_brc_genome_id == Isolate.bv_brc_genome_id)
            .join(Species, Species.species_id == Isolate.species_id)
        ).all()
        for iso, g, species in rows:
            aa = {"gyrA_83": g.gyra83_aa, "gyrA_87": g.gyra87_aa, "parC_80": g.parc80_aa,
                  "parC_84": g.parc84_aa, "recN_346": g.recn346_aa}
            flags = pm.get(iso.bv_brc_genome_id, {})
            res = run_grammar(aa, {"qnrB2": flags.get("qnrB2", 0), "qnrS1": flags.get("qnrS1", 0)}, species, g.recn346_codon)
            if res.prediction == "Indeterminate":
                stats["indeterminate"] += 1
                continue
            iso.septatyper_prediction = SeptaPrediction(res.prediction)
            stats[res.prediction] += 1
        if args.dry_run:
            db.rollback()
        else:
            db.commit()
    print(("DRY RUN " if args.dry_run else "") + f"processed {len(rows)} isolates: {dict(stats)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
