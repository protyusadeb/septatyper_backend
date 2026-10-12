"""Score the SeptaTyper grammar against the lab phenotypes in Table S1 (no BLAST needed).

    python -m scripts.validate_grammar "Table S1.xlsx"

This is concordance on the founding dataset, not independent validation: the grammar
was written with this data in view.
"""
import sys
from pathlib import Path

import pandas as pd

from septatyper.ingest import parse_table_s1
from septatyper.prediction_engine import run_grammar


def main(path: str) -> None:
    df, _ = parse_table_s1(Path(path))
    rows = []
    for r in df.to_dict("records"):
        aa = {"gyrA_83": r["GyrA-83"], "gyrA_87": r["GyrA-87"], "parC_80": r["ParC-80"],
              "parC_84": r["ParC-84"], "recN_346": r["RecN-346"]}
        pm = {"qnrB2": int(r["qnrB2"]), "qnrS1": int(r["qnrS1"])}
        # Table S1 has no codons, so the S. enterica CTG rule cannot be scored here.
        g = run_grammar(aa, pm, r["species"], None)
        rows.append({
            "species": r["species"], "tetratype": "-".join(aa[k] for k in list(aa)[:4]),
            "qnr": int(pm["qnrB2"] or pm["qnrS1"]), "recN": aa["recN_346"],
            "lab": r["lab_phenotype"].value, "pred": g.phenotype,
            "by_h": g.tier1.decided_by == "hydrophobicity",
            "override": g.tier2.override_applied,
        })
    d = pd.DataFrame(rows)
    d["ok"] = d.lab == d.pred
    print(f"Isolates: {len(d)}  concordance: {d.ok.mean():.1%}")
    print(f"Resolved by <H> tie-break: {int(d.by_h.sum())}  Tier 2 overrides applied: {int(d.override.sum())}")
    print(pd.crosstab(d.lab, d.pred))
    print("\nBy species:\n", d.groupby("species").ok.agg(["size", "mean"]).round(3))
    print("\nLargest discordant groups:")
    bad = d[~d.ok].groupby(["species", "tetratype", "qnr", "recN", "lab"]).size().sort_values(ascending=False)
    print(bad.head(8).to_string())


if __name__ == "__main__":
    main(sys.argv[1])
