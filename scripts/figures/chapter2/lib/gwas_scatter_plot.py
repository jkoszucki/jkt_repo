"""
Scatter plot: F1 vs MCC for all pyseer hits (PCI50C50, F1 >= 0.5, MCC >= 0.5).

Every row in pyseer_hits_all.tsv with version == "PCI50C50" and both F1 >= 0.5
and MCC >= 0.5 is plotted as one point.  Axes start at 0.5 (no points below).

Colouring:
  SGNH hydrolase PC  →  style.sgnh_domain_color
                         (any hit in PCI50C50/clusters_functions.tsv with
                          reported_topology_PC == "SGNH hydrolase",
                          tcov >= 0.1, prob >= 0.7)
  All other PCs      →  gray (#bfbfbf)

Labels (only for tail spike / tail fiber PCs):
  For each unique PC find the PHROGS hit with the highest bitscore.
  If its function contains "tail fiber" or "tail spike" (case-insensitive),
  label every point for that PC with the function text.  No other labels are shown.

Output:
  plots_dir/gwas_scatter_f1_mcc.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


_TAIL_FIBER_KW = "tail fiber"
_VERSION       = "PCI50C50"


def _load_annotations(hhsearch_dir: Path) -> pd.DataFrame:
    tsv = hhsearch_dir / _VERSION / "clusters_functions.tsv"
    return pd.read_csv(tsv, sep="\t")


def _topology_pcs(ann: pd.DataFrame, topology: str) -> set[str]:
    """PCs with a given reported_topology_PC (tcov >= 0.1, prob >= 0.7)."""
    mask = (
        (ann["reported_topology_PC"] == topology) &
        (ann["tcov"] >= 0.1) &
        (ann["prob"] >= 0.7)
    )
    return set(ann.loc[mask, "PC"].unique())


_PHROG_OUTLINE_COLORS = {
    "unknown_or_no_hit": ("#cccccc", 0.8),   # pale gray, thin
    "other_phrog":       ("#555555", 2.0),   # gray, thick (same as tail fiber)
}


def _phrog_outline_category(ann: pd.DataFrame, pcs: pd.Index) -> dict[str, str | None]:
    """
    Classify each PC by its best-bitscore PHROG hit (for gray/Other PCs):
      'unknown_function' — function == 'unknown function'
      'no_hit'           — no PHROGS entry at all
      'other_phrog'      — any other known function (tail protein, acetylase, …)
    Returns None for PCs that are tail fiber (handled separately).
    """
    phrogs = ann[ann["db"] == "PHROGS"].copy()
    result: dict[str, str | None] = {}
    for pc in pcs:
        sub = phrogs[phrogs["PC"] == pc]
        if sub.empty:
            result[pc] = "unknown_or_no_hit"
            continue
        best_func = str(sub.loc[sub["bits"].idxmax(), "function"]).lower()
        if _TAIL_FIBER_KW in best_func:
            result[pc] = None  # tail fiber — handled by its own overlay
        elif best_func == "unknown function":
            result[pc] = "unknown_or_no_hit"
        else:
            result[pc] = "other_phrog"
    return result


def _is_tail_fiber(ann: pd.DataFrame, pcs: pd.Index) -> set[str]:
    """Return PCs whose best-bitscore PHROG hit is a tail fiber."""
    phrogs = ann[ann["db"] == "PHROGS"].copy()
    result: set[str] = set()
    for pc in pcs:
        sub = phrogs[phrogs["PC"] == pc]
        if sub.empty:
            continue
        best_func = str(sub.loc[sub["bits"].idxmax(), "function"]).lower()
        if _TAIL_FIBER_KW in best_func:
            result.add(pc)
    return result


def _no_ecod_pcs(ann: pd.DataFrame, pcs) -> set[str]:
    """Return PCs with no ECOD hit at prob >= 0.7."""
    ecod = ann[(ann["db"] == "ECOD") & (ann["prob"] >= 0.7)]
    pcs_with_ecod = set(ecod["PC"].unique())
    return set(pc for pc in pcs if pc not in pcs_with_ecod)


def plot_gwas_scatter(
    pyseer_hits_tsv: Path,
    hhsearch_dir: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Plot F1 vs MCC for all PCI50C50 pyseer hits with F1 >= 0.5 and MCC >= 0.5.

    Args:
        pyseer_hits_tsv:  cfg.input_dir / "gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv"
        hhsearch_dir:     cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/HHSEARCH"
        plots_dir:        scripts/figures/chapter2/plots
        style:            cfg.style
    """
    dpi        = getattr(style, "dpi", 300)
    tick_fs    = getattr(style, "tick_fontsize", 8)
    tick_fw    = getattr(style, "tick_fontweight", "bold")
    label_fs   = getattr(style, "axis_label_fontsize", 12)
    label_fw   = getattr(style, "axis_label_fontweight", "bold")
    sgnh_color   = getattr(style, "sgnh_domain_color", "#c9a227")
    gray_color   = getattr(style, "gray_color", "#bfbfbf")
    pectin_color = "#2ca02c"  # green

    # --- filter pyseer hits ---
    hits = pd.read_csv(pyseer_hits_tsv, sep="\t")
    df = hits[
        (hits["version"] == _VERSION) &
        (hits["F1_score"] >= 0.5) &
        (hits["MCC"] >= 0.5)
    ].copy().reset_index(drop=True)
    print(f"  [scatter] {len(df)} rows ({df['PC'].nunique()} unique PCs, "
          f"{df['locus'].nunique()} loci) after filter.")

    # --- load PCI50C50 annotations ---
    no_ecod_color = "#ffffff"  # white fill for PCs with no ECOD hit

    ann        = _load_annotations(hhsearch_dir)
    sgnh           = _topology_pcs(ann, "SGNH hydrolase")
    pectin         = _topology_pcs(ann, "Pectin lyase-like")
    no_ecod        = _no_ecod_pcs(ann, df["PC"].unique())
    tail_fiber     = _is_tail_fiber(ann, df["PC"].unique())
    phrog_outline  = _phrog_outline_category(ann, df["PC"].unique())

    # SGNH/Pectin take priority; then no ECOD hit → white; otherwise gray
    def _color(pc):
        if pc in sgnh:    return sgnh_color
        if pc in pectin:  return pectin_color
        if pc in no_ecod: return no_ecod_color
        return gray_color

    df["color"]         = df["PC"].map(_color)
    df["tail_fiber"]    = df["PC"].isin(tail_fiber)
    df["phrog_outline"] = df["PC"].map(phrog_outline)

    # --- plot ---
    fig, ax = plt.subplots(figsize=(4.0, 3.8))

    # x=y reference line
    ax.plot([0.47, 1.03], [0.47, 1.03], color="#aaaaaa", lw=0.8, linestyle="--", zorder=1)

    # draw in order: no_ecod (white), gray, pectin, SGNH on top
    zorder_map = {no_ecod_color: 2, gray_color: 3, pectin_color: 4, sgnh_color: 5}
    for color, subset in df.groupby("color", sort=False):
        edge = "#888888" if color == no_ecod_color else "white"
        lw   = 0.6       if color == no_ecod_color else 0.4
        ax.scatter(
            subset["F1_score"], subset["MCC"],
            c=color, s=80, alpha=0.7, zorder=zorder_map.get(color, 2),
            edgecolors=edge, linewidths=lw,
        )

    # legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=sgnh_color,   edgecolor="white",   label="SGNH hydrolase"),
        Patch(facecolor=pectin_color, edgecolor="white",   label="Pectin lyase-like"),
        Patch(facecolor=gray_color,   edgecolor="white",   label="Other ECOD"),
        Patch(facecolor=no_ecod_color, edgecolor="#888888", linewidth=0.6, label="No ECOD hit"),
    ]
    ax.legend(handles=legend_elements, fontsize=tick_fs - 1, frameon=False)

    ax.set_xlim(0.47, 1.03)
    ax.set_ylim(0.47, 1.03)
    label_fs_small = round(label_fs * 0.8)
    ax.set_xlabel("F1 score", fontsize=label_fs_small, fontweight=label_fw)
    ax.set_ylabel("Matthews correlation coefficient (MCC)", fontsize=label_fs_small, fontweight=label_fw)
    ax.tick_params(labelsize=round(tick_fs * 0.8))
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight(tick_fw)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout(pad=1.5)
    out_path = plots_dir / "gwas_scatter_f1_mcc.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  [scatter] → {out_path.name}")


def plot_gwas_phrog_barplot(
    pyseer_hits_tsv: Path,
    hhsearch_dir: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Horizontal bar plot of PHROG function counts for all unique PCs in the
    PCI50C50 filtered set (F1 >= 0.5, MCC >= 0.5), sorted by frequency.

    Args:
        pyseer_hits_tsv:  cfg.input_dir / "gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv"
        hhsearch_dir:     cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/HHSEARCH"
        plots_dir:        scripts/figures/chapter2/plots
        style:            cfg.style
    """
    dpi      = getattr(style, "dpi", 300)
    tick_fs  = getattr(style, "tick_fontsize", 8)
    tick_fw  = getattr(style, "tick_fontweight", "bold")
    label_fs = getattr(style, "axis_label_fontsize", 12)
    label_fw = getattr(style, "axis_label_fontweight", "bold")
    gray_color = getattr(style, "gray_color", "#bfbfbf")

    hits = pd.read_csv(pyseer_hits_tsv, sep="\t")
    pcs = hits[
        (hits["version"] == _VERSION) &
        (hits["F1_score"] >= 0.5) &
        (hits["MCC"] >= 0.5)
    ]["PC"].unique()

    ann    = _load_annotations(hhsearch_dir)
    phrogs = ann[ann["db"] == "PHROGS"].copy()

    func_counts: dict[str, int] = {}
    for pc in pcs:
        sub = phrogs[phrogs["PC"] == pc]
        func = str(sub.loc[sub["bits"].idxmax(), "function"]) if not sub.empty else "no PHROG hit"
        func_counts[func] = func_counts.get(func, 0) + 1

    funcs  = sorted(func_counts, key=func_counts.get, reverse=True)
    counts = [func_counts[f] for f in funcs]

    fig, ax = plt.subplots(figsize=(4.0, max(2.5, len(funcs) * 0.6)))
    ax.barh(funcs, counts, color=gray_color, edgecolor="white", linewidth=0.4)

    ax.set_xlabel("Number of PCs", fontsize=round(label_fs * 0.8), fontweight=label_fw)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.tick_params(axis="y", labelsize=tick_fs)
    ax.tick_params(axis="x", labelsize=round(tick_fs * 0.8))
    for lbl in ax.get_yticklabels():
        lbl.set_fontweight(tick_fw)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    out_path = plots_dir / "gwas_phrog_functions.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  [phrog_bar] → {out_path.name}")
