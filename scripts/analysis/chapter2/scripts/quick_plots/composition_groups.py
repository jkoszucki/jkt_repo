import os
import pandas as pd
import re
from itertools import combinations, product
import matplotlib.pyplot as plt

def parse_set(val):
    if pd.isna(val):
        return frozenset()
    s = str(val).strip()
    tokens = re.split(r"[;,|]\s*|\s+", s)
    toks = [t.strip() for t in tokens if t and t.strip()]
    return frozenset(toks)

def compute_jaccard(s1: frozenset, s2: frozenset) -> float:
    # define Jaccard for empty sets
    if not s1 and not s2:
        return 1.0
    union = s1 | s2
    if len(union) == 0:
        return 0.0
    inter = s1 & s2
    return len(inter) / len(union)

def main(input_csv: str, outdir: str):
    # 1. Create folder
    os.makedirs(outdir, exist_ok=True)

    # 2. Read table & compute composition frequencies
    df = pd.read_csv(input_csv)
    if "monos_unique_all" not in df.columns:
        raise KeyError("'monos_unique_all' column missing in input")

    df["composition_set"] = df["monos_unique_all"].apply(parse_set)

    comp_counts = df["composition_set"].value_counts().sort_values(ascending=False)

    total = comp_counts.sum()
    comp_df = pd.DataFrame({
        "composition": [";".join(sorted(x)) if x else "" for x in comp_counts.index],
        "count": comp_counts.values,
        "frequency": comp_counts.values / total
    })
    comp_df.to_csv(os.path.join(outdir, "composition_freq.csv"), index=False)

    # 3. Pairwise Jaccard in long format
    unique_sets = list(comp_counts.index)
    labels = [";".join(sorted(s)) if s else "" for s in unique_sets]

    records = []
    # Use product to include both directions or combinations for symmetric?
    # If you want symmetric duplicate (A vs B and B vs A), use product.
    # If you want only one direction (A vs B once), use combinations.
    for i, s1 in enumerate(unique_sets):
        for j, s2 in enumerate(unique_sets):
            rec = {
                "composition_1": labels[i],
                "composition_2": labels[j],
                "jaccard": compute_jaccard(s1, s2)
            }
            records.append(rec)

    sim_long_df = pd.DataFrame(records)
    sim_long_df = ytosim_long_df.query('jaccard>=0.8')
    sim_long_df.to_csv(os.path.join(outdir, "composition_sim.csv"), index=False)


    plt.figure(figsize=(5,4))
    sim_long_df["jaccard"].hist(bins=30, edgecolor="black")
    plt.xlabel("Jaccard similarity", fontsize=12, fontweight="bold")
    plt.ylabel("Count", fontsize=12, fontweight="bold")
    plt.title("Distribution of Jaccard similarities", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "composition_jaccard_hist.png"), dpi=300)
    plt.close()

    print("Wrote:")
    print("  ", os.path.join(outdir, "composition_freq.csv"))
    print("  ", os.path.join(outdir, "composition_sim.csv"))


if __name__ == "__main__":
    # Paths
    infile = "/Users/januszkoszucki/Documents/Thesis/chapter2/capsule_structures/analysis/ktypes.csv"
    outdir = "/Users/januszkoszucki/Documents/Thesis/chapter2/capsule_structures/analysis/plots/composition"
    main(infile, outdir)
