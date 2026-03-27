"""
Render phandango-style tree + heatmap PNGs in Python.

For each locus reads the pre-exported files produced by phandango_export:
    plots/phandango/{locus}/_subtree.nwk
    plots/phandango/{locus}/_variants.csv

and renders a cladogram with two colour columns (K-locus presence, best-PC
presence) to plots/phandango/{locus}/{locus}.png.

K-locus column encodes confidence from bacteria_metadata.tsv:
    Perfect / Very high  →  kpam_color         (full blue)
    High                 →  lighter blue
    Good                 →  even lighter blue
    Low                  →  pale red
    Absent               →  #FDFEFE
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
from io import StringIO

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from Bio import Phylo

_ABSENT_COLOR = "#FDFEFE"

# Ordered confidence levels for the legend (highest → lowest)
CONFIDENCE_LEVELS: list[tuple[str, str]] = [
    ("Perfect / Very high", "#1f77b4"),
    ("High",                "#5a9fd4"),
    ("Good",                "#8ec5e8"),
    ("Low",                 "#e8a0a0"),
    ("Absent",              _ABSENT_COLOR),
]

# Internal lookup: metadata value → colour
_CONF_COLOR: dict[str, str] = {
    "Perfect":   "#1f77b4",
    "Very high": "#1f77b4",
    "High":      "#5a9fd4",
    "Good":      "#8ec5e8",
    "Low":       "#e8a0a0",
}


# ---------------------------------------------------------------------------
# Tree helpers
# ---------------------------------------------------------------------------

def _get_leaf_order(root_clade) -> list[str]:
    """DFS traversal — returns leaf names top → bottom (drawing order)."""
    leaves = []
    def _dfs(clade):
        if clade.is_terminal():
            leaves.append(clade.name)
        else:
            for child in clade.clades:
                _dfs(child)
    _dfs(root_clade)
    return leaves


def _compute_positions(root_clade) -> dict[Any, tuple[float, float]]:
    """
    Cladogram layout: all leaves forced to the same x (max depth) so they
    align flush with the heatmap.  Internal nodes sit at their integer depth.
    y = sequential integer for leaves; internal nodes get mean y of children.
    """
    def _max_depth(clade, d):
        if clade.is_terminal():
            return d
        return max(_max_depth(c, d + 1) for c in clade.clades)
    max_d = _max_depth(root_clade, 0)

    pos: dict[Any, tuple[float, float]] = {}
    counter = [0]

    def _assign(clade, depth: int):
        if clade.is_terminal():
            pos[id(clade)] = (max_d, counter[0])
            counter[0] += 1
        else:
            for child in clade.clades:
                _assign(child, depth + 1)
            ys = [pos[id(c)][1] for c in clade.clades]
            pos[id(clade)] = (depth, (min(ys) + max(ys)) / 2.0)

    _assign(root_clade, 0)
    return pos


def _draw_tree(ax: plt.Axes, root_clade, positions: dict, n_leaves: int) -> None:
    """Draw cladogram: root on left, leaves on right, leaf 0 at top."""

    def _draw(clade):
        px, py = positions[id(clade)]
        for child in clade.clades:
            cx, cy = positions[id(child)]
            ax.plot([px, cx], [cy, cy], color="black", lw=0.5)
            _draw(child)
        if not clade.is_terminal():
            child_ys = [positions[id(c)][1] for c in clade.clades]
            ax.plot([px, px], [min(child_ys), max(child_ys)], color="black", lw=0.5)
        else:
            ax.text(px + 0.05, py, clade.name, va="center", ha="left",
                    fontsize=2.3, color="black")

    _draw(root_clade)
    ax.axis("off")
    ax.set_ylim(n_leaves - 0.5, -0.5)
    # Reserve 40% extra x-space after leaf tips so labels don't run into heatmap
    max_d = max(p[0] for p in positions.values())
    ax.set_xlim(0, max_d * 1.4)


# ---------------------------------------------------------------------------
# Heatmap helper
# ---------------------------------------------------------------------------

def _draw_heatmap(
    ax: plt.Axes,
    leaf_order: list[str],
    values: dict[str, float],
    color_present: str,
    color_absent: str,
    conf_colors: dict[str, str] | None = None,
) -> None:
    """
    Draw a single column of filled unit rectangles (one per leaf).

    If *conf_colors* is supplied (genomeID → hex colour), it overrides
    color_present for the K-locus column to encode confidence level.
    Absent cells always use color_absent.
    """
    n = len(leaf_order)
    for i, genome_id in enumerate(leaf_order):
        val = values.get(genome_id, 0)
        if val:
            color = conf_colors.get(genome_id, color_present) if conf_colors else color_present
        elif conf_colors and genome_id in conf_colors:
            # Low-confidence: K_locus matches but GWAS=0
            color = conf_colors[genome_id]
        else:
            color = color_absent
        ax.add_patch(plt.Rectangle((0, n - 1 - i), 1, 1, color=color, linewidth=0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n)
    ax.axis("off")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def render_phandango_png(
    locus: str,
    pc: str,
    plots_dir: Path,
    meta_tsv: Path,
    style=None,
) -> None:
    """
    Render one locus → plots/phandango/{locus}/{locus}.png.

    Args:
        locus:      K-locus name, e.g. "KL1"
        pc:         best-predictor PC column name, e.g. "PC1817"
        plots_dir:  scripts/figures/chapter2/plots
        meta_tsv:   path to bacteria_metadata.tsv
        style:      cfg.style (optional)
    """
    locus_color = getattr(style, "kpam_color",       "#1f77b4")
    pc_color    = getattr(style, "sgnh_domain_color", "#c9a227")

    phandango_dir = plots_dir / "per-locus" / locus
    nwk_path      = phandango_dir / "_subtree.nwk"
    csv_path      = phandango_dir / "_variants.csv"

    if not nwk_path.exists() or not csv_path.exists():
        print(f"  [render] skipping {locus}: files not found in {phandango_dir}")
        return

    # --- parse tree ---
    tree      = Phylo.read(StringIO(nwk_path.read_text()), "newick")
    root      = tree.root
    leaves    = _get_leaf_order(root)
    positions = _compute_positions(root)
    n         = len(leaves)

    # --- load variants ---
    variants = pd.read_csv(csv_path)
    variants["genomeID"] = variants["genomeID"].astype(str)
    locus_vals = dict(zip(variants["genomeID"], variants.get(locus, pd.Series(dtype=float)).fillna(0)))
    pc_vals    = dict(zip(variants["genomeID"], variants.get(pc,    pd.Series(dtype=float)).fillna(0)))

    # --- confidence colour map (genomeID → hex) for K-locus column ---
    conf_colors: dict[str, str] = {}
    if meta_tsv.exists():
        meta = pd.read_csv(
            meta_tsv, sep="\t",
            usecols=["genomeID", "K_locus", "K_locus_confidence"],
        )
        meta["genomeID"] = meta["genomeID"].astype(str)
        locus_meta = meta[meta["K_locus"] == locus]
        for _, row in locus_meta.iterrows():
            color = _CONF_COLOR.get(row["K_locus_confidence"])
            if color:
                conf_colors[row["genomeID"]] = color

    # --- figure layout ---
    # Width: 30% narrower than original (0.7 × 2.4 = 1.68)
    # Columns: tree 70% / each heatmap column 15%  →  width_ratios [7, 1.5, 1.5]
    fig = plt.figure(figsize=(2.2, max(4, n * 0.046)))
    gs  = gridspec.GridSpec(
        2, 3,
        width_ratios=[7, 1.5, 1.5],
        height_ratios=[0.15, 1],
        hspace=0.01, wspace=0.12,
        left=0.01, right=0.99, top=0.97, bottom=0.01,
    )

    ax_h_locus = fig.add_subplot(gs[0, 1])
    ax_h_pc    = fig.add_subplot(gs[0, 2])
    ax_tree    = fig.add_subplot(gs[1, 0])
    ax_locus   = fig.add_subplot(gs[1, 1])
    ax_pc      = fig.add_subplot(gs[1, 2])

    # --- draw ---
    _draw_tree(ax_tree, root, positions, n)
    _draw_heatmap(ax_locus, leaves, locus_vals, locus_color, _ABSENT_COLOR,
                  conf_colors=conf_colors)
    _draw_heatmap(ax_pc, leaves, pc_vals, pc_color, _ABSENT_COLOR)

    for ax_h, label, color in [
        (ax_h_locus, locus, locus_color),
        (ax_h_pc,    pc,    pc_color),
    ]:
        ax_h.set_xlim(0, 1); ax_h.set_ylim(0, 1); ax_h.axis("off")
        ax_h.text(0.5, 0.5, label, ha="center", va="center",
                  fontsize=12, fontweight="bold", color="black",
                  rotation=90, transform=ax_h.transAxes)

    # --- save ---
    dpi = getattr(style, "dpi", 150)
    fig.savefig(phandango_dir / f"{locus}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def render_all_phandango(
    best_predictors_csv: Path,
    plots_dir: Path,
    meta_tsv: Path,
    style=None,
) -> None:
    """
    Render one PNG per row in best_predictors.csv.

    Args:
        best_predictors_csv:  cfg.output_dir / "chapter2" / "best_predictors.csv"
        plots_dir:            scripts/figures/chapter2/plots
        meta_tsv:             path to bacteria_metadata.tsv
        style:                cfg.style
    """
    df = pd.read_csv(best_predictors_csv)
    ok, skipped = 0, []

    for _, row in df.iterrows():
        locus = row["locus"]
        pc    = row["PC"]
        try:
            render_phandango_png(locus, pc, plots_dir, meta_tsv, style)
            ok += 1
        except Exception as exc:
            skipped.append(f"{locus}: {exc}")

    print(f"  Per-locus PNGs rendered: {ok}/{len(df)} → {plots_dir / 'per-locus'}")
    if skipped:
        for s in skipped:
            print(f"    skipped {s}")
