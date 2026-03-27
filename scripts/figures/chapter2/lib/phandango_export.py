"""
Generate per-locus phandango input files (tree + CSV) for the 18 SGNH best predictors.

For each best predictor the source phandango directory already contains a full
_variants.csv with all significant PCs.  This module extracts the two columns
that matter — the K-locus presence/absence and the single best-predictor PC —
plus their :colour columns, and writes lean phandango-ready files to plots_dir.

If the best-predictor PC is absent from the source _variants.csv (it may have
been excluded from the per-locus phandango output), it is read directly from
the mmseqs binary matrix for that clustering version.

Note: Phandango reads colours per column, not per cell, so a single colour is
assigned per column (present vs absent).

Colors from cfg.style:
    locus present  → style.kpam_color      (default #1f77b4)
    PC present     → style.acetylation_color (default #d45515)
    absent         → #FDFEFE

Outputs (one subfolder per locus):
    plots_dir/phandango/{locus}/_subtree.nwk   — copied from source
    plots_dir/phandango/{locus}/_variants.csv  — genomeID + locus + best PC
"""

from __future__ import annotations

import shutil
from io import StringIO

import pandas as pd
from Bio import Phylo
from pathlib import Path

_ABSENT_COLOR = "#FDFEFE"


def _pc_from_binary_matrix(mmseqs_dir: Path, version: str, pc: str) -> pd.Series:
    """Return a Series (index=genomeID) with 0/1 presence for *pc* in *version*."""
    matrix_path = mmseqs_dir / version / "3_binary_matrix.tsv"
    row = pd.read_csv(matrix_path, sep="\t", index_col=0).loc[pc]
    row.name = pc
    return row


def export_phandango(
    best_predictors_csv: Path,
    phandango_src_root: Path,
    mmseqs_dir: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Args:
        best_predictors_csv:  cfg.output_dir / "chapter2" / "best_predictors.csv"
        phandango_src_root:   cfg.input_dir / "gwas/3_GWAS/4_ANALYZE/1_PER_LOCUS"
        mmseqs_dir:           cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/2_MMSEQS"
        plots_dir:            scripts/figures/chapter2/plots
        style:                cfg.style
    """
    locus_present_color = getattr(style, "kpam_color",        "#1f77b4")
    pc_present_color    = getattr(style, "sgnh_domain_color", "#c9a227")

    df = pd.read_csv(best_predictors_csv)
    out_root = plots_dir / "per-locus"
    out_root.mkdir(parents=True, exist_ok=True)

    ok, skipped = 0, []

    for _, row in df.iterrows():
        locus   = row["locus"]
        pc      = row["PC"]
        version = row["version"]
        mode    = row["mode"]

        src_dir      = phandango_src_root / locus / mode / "phandango" / version
        variants_src = src_dir / "_variants.csv"
        tree_src     = src_dir / "_subtree.nwk"

        if not variants_src.exists() or not tree_src.exists():
            skipped.append(locus)
            continue

        variants  = pd.read_csv(variants_src)
        locus_col = locus
        locus_clr = f"{locus}:colour"
        pc_col    = pc
        pc_clr    = f"{pc}:colour"

        if locus_col not in variants.columns:
            skipped.append(f"{locus} (missing locus column)")
            continue

        out = variants[["genomeID", locus_col]].copy()
        out[locus_clr] = out[locus_col].map({1: locus_present_color, 0: _ABSENT_COLOR})

        if pc_col in variants.columns:
            out[pc_col] = variants[pc_col]
        else:
            try:
                pc_series = _pc_from_binary_matrix(mmseqs_dir, version, pc)
            except (KeyError, FileNotFoundError) as e:
                skipped.append(f"{locus} (binary matrix fallback failed: {e})")
                continue
            pc_series = pc_series.reindex(out["genomeID"].values).fillna(0).astype(int)
            out[pc_col] = pc_series.values

        out[pc_clr] = out[pc_col].map({1: pc_present_color, 0: _ABSENT_COLOR})

        # Filter to tree leaves so phandango browser and PNG show identical genomes
        tree_leaves = {
            c.name for c in Phylo.read(StringIO(tree_src.read_text()), "newick").get_terminals()
        }
        out = out[out["genomeID"].isin(tree_leaves)]

        out_dir = out_root / locus
        out_dir.mkdir(exist_ok=True)
        out.to_csv(out_dir / "_variants.csv", index=False)
        shutil.copy(tree_src, out_dir / "_subtree.nwk")
        ok += 1

    print(f"  Per-locus files written for {ok}/{len(df)} loci → {out_root}")
    if skipped:
        print(f"  Skipped: {skipped}")
