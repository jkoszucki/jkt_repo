"""
Per-PC PFAM domain profile with PHROGs tail-fiber/spike annotation.

X-axis: one column per best predictor PC, labelled as:
            {K-locus}
            {PC}({tail fiber}|{tail spike})   if a qualifying hit exists
            {PC}                               otherwise
Y-axis: PFAM domains detected at prob >= 0.90, sorted by frequency.
Cells:  dark gray = present (1), white = absent (0); all cells annotated.
Grid:   minor-tick grid separates every cell.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

PROB_THRESHOLD = 0.90
TAIL_KEYWORDS  = ("tail fiber", "tail spike")


def _short_name(function_str: str) -> str:
    """Return the description part after the first ' | '."""
    if " | " in function_str:
        return function_str.split(" | ", 1)[1].strip()
    return function_str.strip()


def _tail_label(phrogs_df: pd.DataFrame, pc: str) -> str | None:
    """
    Return the first matching tail keyword for this PC, or None.
    Returns 'tail fiber' or 'tail spike'.
    """
    hits = phrogs_df.loc[phrogs_df["PC"] == pc, "function"].str.lower()
    for kw in TAIL_KEYWORDS:
        if hits.str.contains(kw).any():
            return kw
    return None


def _x_label(locus: str, pc: str, tail: str | None) -> str:
    if tail:
        return f"{locus}\n{pc}\n({tail})"
    return f"{locus}\n{pc}"


def plot_functions_heatmap(
    functions_tsv: Path,
    best_predictors_csv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Write functions_heatmap.png.

    Args:
        functions_tsv:       path to best_predictors_functions.tsv
        best_predictors_csv: path to best_predictors.csv (for PC → locus mapping)
        plots_dir:           output directory for the plot
        style:               cfg.style
    """
    dpi      = getattr(style, "dpi", 300)
    tick_fs  = getattr(style, "tick_fontsize", 8)
    tick_fw  = getattr(style, "tick_fontweight", "bold")
    label_fs = getattr(style, "axis_label_fontsize", 12)
    label_fw = getattr(style, "axis_label_fontweight", "bold")

    # --- Load & filter ---
    df   = pd.read_csv(functions_tsv, sep="\t")
    df09 = df[df["prob"] >= PROB_THRESHOLD].copy()

    locus_map = (
        pd.read_csv(best_predictors_csv)[["PC", "locus"]]
        .set_index("PC")["locus"]
        .to_dict()
    )

    pfam   = df09[df09["db"] == "PFAM"][["PC", "function"]].copy()
    phrogs = df09[df09["db"] == "PHROGS"][["PC", "function"]].copy()

    pfam["function"] = pfam["function"].map(_short_name)
    pfam = pfam.drop_duplicates(subset=["PC", "function"])

    # --- Build presence/absence matrix ---
    pfam["present"] = 1
    matrix = (
        pfam.pivot_table(index="function", columns="PC",
                         values="present", aggfunc="max")
        .fillna(0).astype(int)
    )

    # Sort rows by frequency descending
    row_order = matrix.sum(axis=1).sort_values(ascending=False).index.tolist()

    # Sort columns by Jaccard similarity of PFAM profiles (hierarchical clustering)
    # matrix.T: rows = PCs, cols = PFAM domains — Jaccard distance on binary vectors
    dist = pdist(matrix.values.T, metric="jaccard")
    linkage_matrix = linkage(dist, method="average")
    col_idx = leaves_list(linkage_matrix)
    col_order = [matrix.columns[i] for i in col_idx]

    # Two manual groups, each sorted by K-locus number
    def _locus_num(pc):
        return int("".join(filter(str.isdigit, locus_map.get(pc, "0"))) or "0")

    RIGHT_LOCI = {"KL2", "KL7", "KL24", "KL30", "KL48", "KL57"}
    left_pcs  = sorted(
        [pc for pc in col_order if locus_map.get(pc) not in RIGHT_LOCI],
        key=_locus_num,
    )
    right_pcs = sorted(
        [pc for pc in col_order if locus_map.get(pc) in RIGHT_LOCI],
        key=_locus_num,
    )
    col_order = left_pcs + right_pcs

    matrix = matrix.reindex(index=row_order, columns=col_order)

    # Build X-axis labels
    x_labels = [
        _x_label(locus_map.get(pc, pc), pc, _tail_label(phrogs, pc))
        for pc in col_order
    ]

    # --- Plot ---
    n_rows, n_cols = matrix.shape
    cell_w, cell_h = 0.70, 0.42
    fig_w = max(6, n_cols * cell_w + 4)
    fig_h = max(4, n_rows * cell_h + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    gray_cmap = mcolors.LinearSegmentedColormap.from_list(
        "binary_gray", ["#ffffff", "#2b2b2b"]
    )
    ax.imshow(matrix.values, aspect="auto", cmap=gray_cmap, vmin=0, vmax=1)

    # Major ticks for labels
    ax.set_xticks(range(n_cols))
    ax.set_yticks(range(n_rows))
    ax.set_xticklabels(x_labels, rotation=45, ha="right",
                       fontsize=tick_fs, fontweight=tick_fw)
    ax.set_yticklabels(row_order, fontsize=tick_fs, fontweight=tick_fw)

    # Minor ticks for grid lines between cells
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="#aaaaaa", linewidth=0.6)
    ax.tick_params(which="minor", length=0)

    ax.set_xlabel("Best predictor (K-locus)", labelpad=8,
                  fontsize=label_fs, fontweight=label_fw)
    ax.set_ylabel("PFAM domain", labelpad=8,
                  fontsize=label_fs, fontweight=label_fw)

    plt.tight_layout()
    out_path = plots_dir / "functions_heatmap.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out_path}")
