"""
Step 5: Add enzyme_tmscore_05 column to pyseer_hits_sslbh_pvalcor005.tsv.

Reads no-ecod-reported-topology-tm-align.tsv, filters rows where
tm_score_query >= 0.5 OR tm_score_ref >= 0.5, and adds a boolean column
enzyme_tmscore_05 to the pyseer table (True for matching clustering_level/PC pairs).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def merge_tmalign_into_pyseer(tmalign_tsv: Path, acetyl_dir: Path) -> None:
    """
    Add enzyme_tmscore_05 (bool) column to pyseer_hits_sslbh_pvalcor005.tsv.

    Args:
        tmalign_tsv: acetyl-gwas/no-ecod-reported-topology-tm-align.tsv
        acetyl_dir:  output_dir/acetyltransferase/
    """
    pyseer_tsv = acetyl_dir / "acetyl-gwas" / "pyseer_hits_sslbh_pvalcor005.tsv"

    if not tmalign_tsv.exists():
        print("  [skip] TM-align results not found — skipping")
        return

    tm = pd.read_csv(tmalign_tsv, sep="\t")
    if tm.empty:
        print("  [skip] TM-align results empty — skipping")
        return

    positive = tm[(tm["tm_score_query"] >= 0.5) | (tm["tm_score_ref"] >= 0.5)]
    positive_pairs = set(zip(positive["clustering_level"], positive["pc"]))
    print(f"  TM-score ≥ 0.5 pairs: {len(positive_pairs)}")

    pyseer = pd.read_csv(pyseer_tsv, sep="\t")

    # Drop any previously added TM-align columns (idempotent re-run)
    drop = [c for c in pyseer.columns
            if c.startswith("PROTEIN0") or c == "enzyme_tmscore_05"]
    pyseer = pyseer.drop(columns=drop, errors="ignore")

    pyseer["enzyme_tmscore_05"] = [
        (row["clustering_level"], row["PC"]) in positive_pairs
        for _, row in pyseer.iterrows()
    ]

    pyseer.to_csv(pyseer_tsv, sep="\t", index=False)
    n = pyseer["enzyme_tmscore_05"].sum()
    print(f"  enzyme_tmscore_05=True: {n} rows")
    print(f"  → {pyseer_tsv}")
