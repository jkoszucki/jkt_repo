#!/usr/bin/env python3
# Scatter: per clustering level, X=LQ best, Y=HQ best (one point per K locus).
# Labels: annotate outliers (per panel) by K-locus.

import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Paths
results_dir = Path("../results")
in_path     = results_dir / "bestf1.tsv"
out_png     = results_dir / "bestf1_scatter.png"
out_pdf     = results_dir / "bestf1_scatter.pdf"

# Grid definition: columns = coverage (50, 80); rows = identity (00, 50, 80)
id_order  = [0, 50, 80]   # rows (top→bottom)
cov_order = [50, 80]      # columns (left→right)
required_versions = ["PCI00C50", "PCI00C80",
                     "PCI50C50", "PCI50C80",
                     "PCI80C50", "PCI80C80"]

# ---------- Load ----------
df = pd.read_csv(in_path, sep="\t", dtype=str)

# Required columns
for col in ["version", "locus", "quality"]:
    if col not in df.columns:
        raise ValueError(f"Missing required column: {col}")

# Numeric F1 helper (keep original if present)
if "F1_score_num" not in df.columns:
    if "F1_score" not in df.columns:
        raise ValueError("Need F1_score or F1_score_num column in the input.")
    df["F1_score_num"] = pd.to_numeric(df["F1_score"], errors="coerce")
else:
    df["F1_score_num"] = pd.to_numeric(df["F1_score_num"], errors="coerce")

# Keep only HQ/LQ
df = df[df["quality"].isin(["hq", "lq"])].copy()

# ---------- Helpers ----------
def parse_ic(ver: str):
    # "PCI50C80" -> (identity=50, coverage=80)
    m = re.match(r"^PCI(\d{2})C(\d{2})$", str(ver))
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))

def grid_pos(ver: str):
    # row by identity, column by coverage
    ident, cov = parse_ic(ver)
    if ident is None:
        return None
    try:
        r = id_order.index(ident)   # row = identity
        c = cov_order.index(cov)    # col = coverage
    except ValueError:
        return None
    return (r, c)

# ---------- Figure ----------
fig, axes = plt.subplots(3, 2, figsize=(8.5, 10), sharex=True, sharey=True)

# Global axis limits for scores in [0,1]
xmin = ymin = 0.0
xmax = ymax = 1.0

# ---------- Plot each panel ----------
for ver in required_versions:
    pos = grid_pos(ver)
    if pos is None:
        continue
    r_idx, c_idx = pos
    ax = axes[r_idx, c_idx]

    sub = df[df["version"] == ver].copy()

    # Pivot to get per-locus columns 'lq' and 'hq'
    pvt = (
        sub.pivot_table(
            index="locus",
            columns="quality",
            values="F1_score_num",
            aggfunc="max"
        )
        .reindex(columns=["lq", "hq"])
    )
    # Keep loci that have BOTH values
    pvt = pvt.dropna(subset=["lq", "hq"], how="any")

    x = pvt["lq"].astype(float).values
    y = pvt["hq"].astype(float).values
    loci = pvt.index.to_numpy()

    # Scatter
    ax.scatter(x, y, s=30, alpha=0.9)

    # Diagonal reference (HQ=LQ)
    ax.plot([xmin, xmax], [ymin, ymax], linestyle="--", linewidth=1)

    # Limits & aspect
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal', adjustable='box')

    # -------- Pearson r annotation (replaces "n=..." panel label) --------
    if len(x) >= 2 and np.std(x) > 0 and np.std(y) > 0:
        r_val = np.corrcoef(x, y)[0, 1]
        label = f"r = {r_val:.2f}"
    else:
        label = "r = NA"

    ax.text(
        0.02, 0.98, label,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=9,
    )
    # --------------------------------------------------------------------

    # -------- Outlier labeling (IQR rule on Δ = HQ - LQ) --------
    if len(x) > 0:
        delta = y - x
        q1, q3 = np.percentile(delta, [25, 75])
        iqr = q3 - q1
        if iqr == 0:
            # Fallback: label top 3 by absolute deviation from diagonal
            k = min(3, len(delta))
            idx = np.argsort(np.abs(delta))[-k:]
            mask = np.zeros(len(delta), dtype=bool)
            mask[idx] = True
        else:
            lower = q1 - 1.5 * iqr
            upper = q3 + 1.5 * iqr
            mask = (delta < lower) | (delta > upper)

        # Annotate selected points with K-locus labels
        for xi, yi, lab, is_out in zip(x, y, loci, mask):
            if is_out:
                ax.annotate(
                    str(lab),
                    xy=(xi, yi),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=7.5,
                    fontweight="bold",
                )

# Shared labels and headers
fig.supxlabel(
    "Top F1 score\nHigh and low quality prophages (default)",
    fontsize=12,
    fontweight="bold",
    ha="center",
    va="center",
    y=0.065,
)
fig.supylabel(
    "Top F1 score\nHigh quality prophages",
    fontsize=12,
    fontweight="bold",
    ha="center",
    va="center",
    x=0.1,
)

# Column headers (coverage) and row headers (identity on the RIGHT)
for c2, cov in enumerate(cov_order):
    axes[0, c2].annotate(
        f"Coverage: {cov}%",
        xy=(0.5, 1.05),
        xycoords="axes fraction",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )

for r2, ident in enumerate(id_order):
    axes[r2, -1].annotate(
        f"Identity: {ident}%",
        xy=(1.02, 0.5),
        xycoords="axes fraction",
        ha="left",
        va="center",
        fontsize=10,
        fontweight="bold",
    )

fig.tight_layout(rect=[0.08, 0.05, 0.95, 0.95])
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
