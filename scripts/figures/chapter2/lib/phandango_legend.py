"""Generate a PNG legend for the phandango colour scheme."""

from __future__ import annotations

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from pathlib import Path

from phandango_export import _ABSENT_COLOR
from phandango_render import CONFIDENCE_LEVELS


def plot_phandango_legend(plots_dir: Path, style=None) -> None:
    pc_present_color = getattr(style, "sgnh_domain_color", "#c9a227")
    dpi              = getattr(style, "dpi", 300)

    sections = [
        ("K-locus confidence", CONFIDENCE_LEVELS),
        ("SGNH predictor",     [("Present", pc_present_color), ("Absent", _ABSENT_COLOR)]),
    ]

    # Flatten into a single list: section header rows + entry rows
    rows: list[tuple[str, str | None]] = []   # (label, color_or_None)
    for title, entries in sections:
        rows.append((title, None))             # section header
        for label, color in entries:
            rows.append((label, color))

    row_h   = 0.30                             # inches per row
    fig_h   = len(rows) * row_h + 0.1
    fig, ax = plt.subplots(figsize=(2.2, fig_h))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(rows))

    for i, (label, color) in enumerate(rows):
        y = len(rows) - 1 - i                 # top → bottom
        if color is None:
            # Section header
            ax.text(0.02, y + 0.5, label, fontsize=7, fontweight="bold",
                    va="center", ha="left")
        else:
            patch = mpatches.FancyBboxPatch(
                (0.02, y + 0.1), 0.18, 0.65,
                boxstyle="round,pad=0.01",
                facecolor=color,
                edgecolor="0.4" if color == _ABSENT_COLOR else "none",
                linewidth=0.5,
            )
            ax.add_patch(patch)
            ax.text(0.26, y + 0.5, label, fontsize=7, va="center", ha="left")

    fig.tight_layout(pad=0.2)
    out = plots_dir / "per-locus_legend.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    pdf_out = plots_dir / "per-locus_legend.pdf"
    fig.savefig(pdf_out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}, {pdf_out.name}")
