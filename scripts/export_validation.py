"""Score the grammar on Table S1 and write the JSON the Module D dashboard displays.

    python -m scripts.export_validation "Table S1.xlsx" ../septatyper_frontend/data/validation.json

The grammar was developed with this dataset in view, so this is concordance on the founding data,
not independent validation. The S. enterica CTG rule needs codons, which Table S1 lacks, so it is not scored.
"""
import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd

from septatyper.ingest import parse_table_s1
from septatyper.prediction_engine import run_grammar


def main(xlsx: str, out: str) -> None:
    df, _ = parse_table_s1(Path(xlsx))
    rows = []
    for r in df.to_dict("records"):
        aa = {"gyrA_83": r["GyrA-83"], "gyrA_87": r["GyrA-87"], "parC_80": r["ParC-80"],
              "parC_84": r["ParC-84"], "recN_346": r["RecN-346"]}
        g = run_grammar(aa, {"qnrB2": int(r["qnrB2"]), "qnrS1": int(r["qnrS1"])}, r["species"], None)
        rules = g.tier2.rules if g.tier2 else []
        rows.append({
            "species": r["species"], "lab": r["lab_phenotype"].value, "pred": g.phenotype,
            "qnr_rule": "qnr_present" in rules, "baseline": g.tier1.baseline if g.tier1 else None,
        })
    d = pd.DataFrame(rows)
    d["ok"] = d.lab == d.pred
    tp = int(((d.lab == "Resistant") & (d.pred == "Resistant")).sum())
    tn = int(((d.lab == "Susceptible") & (d.pred == "Susceptible")).sum())
    fp = int(((d.lab == "Susceptible") & (d.pred == "Resistant")).sum())
    fn = int(((d.lab == "Resistant") & (d.pred == "Susceptible")).sum())
    # isolates the qnr override moved from a Susceptible baseline to Resistant
    qnr = d[d.qnr_rule]
    by_species = [
        {"species": sp, "n": int(len(x)), "concordance": round(float(x.ok.mean()), 4)}
        for sp, x in d.groupby("species")
    ]
    result = {
        "source": "Table S1, 2,734 isolates with laboratory ciprofloxacin phenotype",
        "generated": date.today().isoformat(),
        "n_isolates": int(len(d)),
        "concordance": round(float(d.ok.mean()), 4),
        "sensitivity": round(tp / (tp + fn), 4),
        "specificity": round(tn / (tn + fp), 4),
        "confusion": {"tp": tp, "tn": tn, "fp": fp, "fn": fn},
        "by_species": by_species,
        "qnr_override": {
            "applied": int(len(qnr)),
            "lab_resistant": int((qnr.lab == "Resistant").sum()),
            "lab_susceptible": int((qnr.lab == "Susceptible").sum()),
            "lab_susceptible_by_species": {sp: int(((x.lab == "Susceptible")).sum()) for sp, x in qnr.groupby("species")},
        },
        "not_scored": "S. enterica F-D-S-E with RecN codon ctg (Table S1 has no codon data)",
    }
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
