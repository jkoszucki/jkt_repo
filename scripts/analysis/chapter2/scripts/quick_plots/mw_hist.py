import os
import pandas as pd
import re
from itertools import combinations, product
import matplotlib.pyplot as plt

input_csv = '/Users/januszkoszucki/Documents/Thesis/chapter2/capsule_structures/analysis/ktypes.csv'
output_path = '/Users/januszkoszucki/Documents/Thesis/chapter2/capsule_structures/analysis/plots/mw_hist.png'

df = pd.read_csv(input_csv)

series = pd.to_numeric(df['mw_struct_kda'], errors="coerce").dropna()

fig, ax = plt.subplots(figsize=(4, 3))
ax.hist(series, bins=40, edgecolor="black")
ax.set_xlabel("mw_struct_kda", fontsize=14, fontweight="bold")
ax.set_ylabel("count", fontsize=14, fontweight="bold")

# Optionally make tick labels bold too:
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(8)
    label.set_fontweight("bold")


fig.savefig(output_path, bbox_inches="tight", dpi=200)
