import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import PercentFormatter


input_csv = '/Users/januszkoszucki/Documents/Thesis/chapter2/capsule_structures/analysis/ktypes_sim.csv'
output_path = '/Users/januszkoszucki/Documents/Thesis/chapter2/capsule_structures/analysis/plots/branch_core_heatmap.png'

# Load
sim_long_df = pd.read_csv(input_csv)

# Keep rows with both values present
sim_long_df = sim_long_df.dropna(subset=["weighted_core", "weighted_total"])
n_total = len(sim_long_df)
if n_total == 0:
    raise ValueError("No data points after dropping NaNs in 'weighted_core' and 'weighted_total'.")

plt.figure(figsize=(5, 4))

# First draw hexbin to get counts per bin
hb = plt.hexbin(
    sim_long_df["weighted_core"].values,
    sim_long_df["weighted_total"].values,
    gridsize=15,
    cmap="Reds",
    mincnt=1  # only show bins that have at least one point
)

# Convert counts to fractions of all comparisons
counts = hb.get_array()
fractions = counts / float(n_total)

# Update the hexbin colors to reflect fractions (0..1)
hb.set_array(fractions)
hb.set_norm(Normalize(vmin=0.0, vmax=fractions.max() if fractions.size else 1.0))

# Percent colorbar
cbar = plt.colorbar(hb, format=PercentFormatter(xmax=1))
cbar.set_label("\nFraction of all comparisons", fontsize=10, fontweight="bold")

plt.xlabel("\nsimilarity score (core)", fontsize=10, fontweight="bold")
plt.ylabel("\nsimilarity score (core + branch)", fontsize=10, fontweight="bold")
plt.title(f"", fontsize=14, fontweight="bold")

# Styling
plt.tick_params(axis="both", which="major", labelsize=12, width=1.5)
plt.axvline(x=0.9, color="gray", linestyle="--", linewidth=1.5)

ax = plt.gca()
for label in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
    label.set_fontweight("bold")
    label.set_fontsize(8)

# plt.axline((0, 0), slope=1, color='black', alpha=0.5, linestyle='-', linewidth=0.5)


plt.tight_layout()
plt.savefig(output_path, dpi=300)
plt.close()
