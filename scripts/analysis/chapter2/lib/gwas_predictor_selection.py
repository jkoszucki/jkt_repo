"""
Select best GWAS predictors per K-locus and annotate with representative sequences.

Pipeline:
1. Load pyseer_hits_all.tsv
2. Filter: precision >= 0.8
3. For each locus, retain the single (PC, version) with highest F1 * MCC
4. Filter: F1 >= 0.5 and MCC >= 0.5
5. Map representative sequence (first sequence of alignment FASTA) per (PC, version)
6. Write best_predictors.csv
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
                    break  # second record — stop
                protein_id = line[1:]  # strip leading ">"
                found = True
            elif found:
                seq_lines.append(line)
    return protein_id, "".join(seq_lines)


def select_best_predictors(table_path: Path, mmseqs_dir: Path, output_dir: Path) -> pd.DataFrame:
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {table_path} ...")
    df = pd.read_csv(table_path, sep="\t")
    print(f"  {len(df):,} rows loaded.")

    # Step 1 — filter precision >= 0.8
    df = df[df["precision"] >= 0.8].copy()
    print(f"  {len(df):,} rows after precision >= 0.8.")

    # Step 2 — score each row, pick best (PC, version) per locus
    df["F1_MCC"] = df["F1_score"] * df["MCC"]
    best = df.loc[df.groupby("locus")["F1_MCC"].idxmax()].copy()
    print(f"  {len(best):,} loci after selecting best (PC, version) per locus.")

    # Step 3 — filter F1 >= 0.5 and MCC >= 0.5
    best = best[(best["F1_score"] >= 0.5) & (best["MCC"] >= 0.5)].copy()
    print(f"  {len(best):,} loci after F1 >= 0.5 and MCC >= 0.5.")

    best = best.sort_values("locus").reset_index(drop=True)

    # Step 4 — map representative sequences
    print("Mapping representative sequences ...")
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

    return best
