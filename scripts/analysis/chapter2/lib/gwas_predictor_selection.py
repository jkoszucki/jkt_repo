"""
Select best GWAS predictors per K-locus and annotate with representative sequences.

Pipeline:
1. Load HHSEARCH annotation files; identify (version, PC) with tcov >= 0.10
   and reported_topology_PC == 'SGNH hydrolase'
2. Load pyseer_hits_all.tsv; filter precision >= 0.6
3. Inner-join with SGNH-qualified (version, PC) pairs
4. For each locus, retain the single (PC, version) with highest F1_score
5. Map representative sequence (first sequence of alignment FASTA) per (PC, version)
6. Write best_predictors.csv
7. Collect all hhsearch rows for selected (version, PC) pairs; write best_predictors_functions.tsv
"""

import pandas as pd
from pathlib import Path


def _read_first_fasta(fasta_path: Path) -> tuple[str, str]:
    """Return (protein_id, sequence) for the first record in a FASTA file."""
    protein_id = ""
    seq_lines = []
    found = False
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if found:
                    break
                protein_id = line[1:]
                found = True
            elif found:
                seq_lines.append(line)
    return protein_id, "".join(seq_lines)


def _load_sgnh_pcs(hhsearch_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load all clustering-level annotation files.

    Returns:
        sgnh_pcs:  (version, PC) pairs with tcov >= 0.10 and SGNH hydrolase fold
        all_hits:  full annotation table across all clustering levels
    """
    dfs = []
    for version_dir in sorted(hhsearch_dir.iterdir()):
        tsv = version_dir / "clusters_functions.tsv"
        if not tsv.exists():
            continue
        df = pd.read_csv(tsv, sep="\t")
        dfs.append(df)

    all_hits = pd.concat(dfs, ignore_index=True)
    print(f"  {len(all_hits):,} annotation rows loaded across {len(dfs)} clustering levels.")

    sgnh = all_hits[
        (all_hits["tcov"] >= 0.10) &
        (all_hits["reported_topology_PC"] == "SGNH hydrolase")
    ][["version", "PC"]].drop_duplicates()
    print(f"  {len(sgnh):,} unique (version, PC) pairs with SGNH hydrolase and tcov >= 10%.")
    return sgnh, all_hits


def select_best_predictors(
    table_path: Path,
    mmseqs_dir: Path,
    hhsearch_dir: Path,
    output_dir: Path,
) -> pd.DataFrame:
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1 — identify SGNH-qualified (version, PC) pairs
    print("Step 1: Loading SGNH hydrolase annotations ...")
    sgnh_pcs, all_hits = _load_sgnh_pcs(hhsearch_dir)

    # Step 2 — load GWAS hits, filter precision >= 0.6
    print(f"Step 2: Loading {table_path.name} ...")
    df = pd.read_csv(table_path, sep="\t")
    print(f"  {len(df):,} rows loaded.")
    df = df[df["precision"] >= 0.6].copy()
    print(f"  {len(df):,} rows after precision >= 0.6.")

    # Step 3 — keep only SGNH-annotated PCs
    df = df.merge(sgnh_pcs, on=["version", "PC"], how="inner")
    print(f"  {len(df):,} rows after SGNH hydrolase filter.")

    # Step 4 — pick highest F1_score per locus (across all clustering levels)
    best = df.loc[df.groupby("locus")["F1_score"].idxmax()].copy()
    best = best.sort_values("locus").reset_index(drop=True)
    print(f"  {len(best):,} loci after selecting best F1 per locus.")

    # Step 5 — map representative sequences
    print("Step 5: Mapping representative sequences ...")
    protein_ids, sequences, missing = [], [], 0
    for _, row in best.iterrows():
        fasta_path = mmseqs_dir / row["version"] / "alignments" / f"{row['PC']}.fasta"
        if not fasta_path.exists():
            print(f"  WARNING: {fasta_path} not found")
            protein_ids.append("")
            sequences.append("")
            missing += 1
            continue
        pid, seq = _read_first_fasta(fasta_path)
        protein_ids.append(pid)
        sequences.append(seq)
    print(f"  Done. Missing: {missing}.")

    best["representative_protein_id"] = protein_ids
    best["representative_sequence"]   = sequences

    out_path = output_dir / "best_predictors.csv"
    best.to_csv(out_path, index=False)
    print(f"Saved → {out_path}")

    # Step 7 — dump all hhsearch rows for the selected (version, PC) pairs
    #           drop ALANDB; keep only hits with prob >= 0.70
    print("Step 7: Exporting hhsearch annotations for best predictors ...")
    selector = best[["version", "PC"]].drop_duplicates()
    best_functions = all_hits.merge(selector, on=["version", "PC"], how="inner")
    best_functions = best_functions[
        (best_functions["db"] != "ALANDB") &
        (best_functions["prob"] >= 0.70)
    ]
    func_path = output_dir / "best_predictors_functions.tsv"
    best_functions.to_csv(func_path, sep="\t", index=False)
    print(f"  {len(best_functions):,} rows → {func_path}")

    return best
