"""
Build all-vs-all blastp similarity network for SGNH GWAS hits + experimental proteins.

Outputs:
    sgnh_network_nodes.tsv  — one row per protein, used as Cytoscape node table
    sgnh_network_edges.tsv  — pairwise BLAST edges with coverage interval labels
"""

import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO


# ---------------------------------------------------------------------------
# Experimental protein metadata
# Substrate/specificity and category for each known experimental protein.
# ---------------------------------------------------------------------------

_EXPERIMENTAL_META = {
    "164_08_KL2":     {"specificity": "4-nitrophenyl acetate", "category": "EXPERIMENTAL_4NPA"},
    "174_38_KL55":    {"specificity": "4-nitrophenyl acetate", "category": "EXPERIMENTAL_4NPA"},
    "YP_004782195.1": {"specificity": "acetylated O-antigen LPS", "category": "EXPERIMENTAL_LPS"},
    "QNO11465.1_K26": {"specificity": "KL26",                  "category": "EXPERIMENTAL_KL26"},
}


# ---------------------------------------------------------------------------
# Node assembly
# ---------------------------------------------------------------------------

def _load_gwas_nodes(best_predictors_csv: Path) -> pd.DataFrame:
    """Load SGNH best-hit proteins from best_predictors.csv."""
    df = pd.read_csv(best_predictors_csv)
    df = df[df["representative_sequence"].notna() & (df["representative_sequence"] != "")]
    nodes = pd.DataFrame({
        "proteinID":  df["representative_protein_id"].values,
        "source":     "PREDICTION",
        "specificity": df["locus"].values,
        "seq":        df["representative_sequence"].values,
    })
    nodes["label1"] = nodes["proteinID"] + "_" + nodes["specificity"]
    nodes["label2"] = nodes["specificity"]
    return nodes


def _load_experimental_nodes(fasta_path: Path) -> pd.DataFrame:
    """Load experimental protein(s) from a FASTA file."""
    rows = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        meta = _EXPERIMENTAL_META.get(rec.id, {})
        specificity = meta.get("specificity")
        category    = meta.get("category", "EXPERIMENTAL")
        rows.append({
            "proteinID":   rec.id,
            "source":      category,
            "specificity": specificity,
            "seq":         str(rec.seq),
            "label1":      rec.id,
            "label2":      specificity if specificity else rec.id,
        })
    return pd.DataFrame(rows)


def build_nodes(best_predictors_csv: Path, experimental_fastas: list[Path]) -> pd.DataFrame:
    parts = [_load_gwas_nodes(best_predictors_csv)]
    for fasta in experimental_fastas:
        parts.append(_load_experimental_nodes(fasta))
    nodes = pd.concat(parts, ignore_index=True)
    nodes.index = nodes.index + 1
    nodes["query"] = nodes.index.astype(str) + "_" + nodes["label1"]
    return nodes[["query", "proteinID", "source", "specificity", "label1", "label2", "seq"]]


# ---------------------------------------------------------------------------
# Edge generation (adapted from FIGURE4_PANELB / generate_edges)
# ---------------------------------------------------------------------------

def _sanitize_id(s: str) -> str:
    """Replace dots with _DOT_ so blastp doesn't truncate IDs."""
    return s.replace(".", "_DOT_")


def _restore_id(s: str) -> str:
    return s.replace("_DOT_", ".")


def _write_fasta(nodes: pd.DataFrame, fasta_path: Path) -> None:
    """Write FASTA using sanitized blast_id (no dots)."""
    with open(fasta_path, "w") as f:
        for _, row in nodes.iterrows():
            f.write(f">{row['blast_id']}\n{row['seq']}\n")


def generate_edges(
    nodes: pd.DataFrame,
    tmp_dir: Path,
    min_cov: float = 0,
    min_pident: float = 0,
    max_evalue: float = 1e-3,
    intervals: list = None,   # [low_bound, mid, high] for coverage %
) -> pd.DataFrame:
    """
    Run all-vs-all blastp on nodes['query'] / nodes['seq'], filter by coverage
    and percent identity, and assign a coverage interval label to each edge.

    intervals default: [0, 50, 80]  → low (0–50 %), medium (50–80 %), high (80–100 %)
    Dots in query IDs are temporarily replaced with _DOT_ for BLAST compatibility.
    """
    if intervals is None:
        intervals = [0, 50, 80]
    bottom, mid, top = intervals

    tmp_dir.mkdir(parents=True, exist_ok=True)
    fasta_path  = tmp_dir / "sequences.fasta"
    blast_raw   = tmp_dir / "blastp_raw.tsv"

    # sanitize IDs for BLAST, keep mapping back to original query IDs
    nodes = nodes.copy()
    nodes["blast_id"] = nodes["query"].apply(_sanitize_id)
    blast_to_query = dict(zip(nodes["blast_id"], nodes["query"]))

    _write_fasta(nodes, fasta_path)

    subprocess.run(
        ["makeblastdb", "-in", str(fasta_path), "-dbtype", "prot"],
        check=True, capture_output=True,
    )
    subprocess.run(
        ["blastp",
         "-query", str(fasta_path), "-db", str(fasta_path),
         "-out", str(blast_raw),
         "-outfmt", "6 qseqid sseqid evalue pident bitscore qlen qstart qend slen sstart send"],
        check=True, capture_output=True,
    )

    cols = ["query", "target", "evalue", "pident", "bitscore",
            "qlen", "qstart", "qend", "slen", "sstart", "send"]
    df = pd.read_csv(blast_raw, sep="\t", header=None, names=cols)

    # restore original IDs (undo dot sanitization)
    df["query"]  = df["query"].map(blast_to_query).fillna(df["query"].apply(_restore_id))
    df["target"] = df["target"].map(blast_to_query).fillna(df["target"].apply(_restore_id))

    df["qcov"] = np.round((df["qend"] - df["qstart"] + 1) / df["qlen"], 2)
    df["scov"] = np.round((df["send"] - df["sstart"] + 1) / df["slen"], 2)

    # deduplicate symmetric hits
    df["pair"] = df.apply(lambda r: "-".join(sorted([r["query"], r["target"]])), axis=1)
    df = df.drop_duplicates(subset="pair")

    # filter
    is_self      = df["query"] == df["target"]
    filt_cov     = (df["qcov"] >= min_cov) & (df["scov"] >= min_cov)
    filt_pident  = df["pident"] >= min_pident
    filt_eval    = df["evalue"] <= max_evalue
    df = df.loc[filt_cov & filt_pident & filt_eval & ~is_self].copy()

    # coverage interval
    cov = df[["qcov", "scov"]].min(axis=1) * 100
    df["coverage_interval"] = pd.cut(
        cov,
        bins=[bottom, mid, top, 100],
        labels=[f"low ({bottom}–{mid}%)", f"medium ({mid}–{top}%)", f"high ({top}–100%)"],
        right=True, include_lowest=True,
    )

    # add singletons
    blast_nodes = set(df["query"]) | set(df["target"])
    singletons = [n for n in nodes["query"] if n not in blast_nodes]
    if singletons:
        sing_df = pd.DataFrame({"query": singletons, "target": singletons,
                                "coverage_interval": f"high ({top}–100%)"})
        df = pd.concat([df, sing_df], ignore_index=True)

    df["interaction"] = "A"
    return df


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def compute_sgnh_network(
    best_predictors_csv: Path,
    experimental_fastas: list[Path],
    output_dir: Path,
    tmp_dir: Path,
) -> None:
    """
    Build nodes + edges for the SGNH similarity network and write TSV files.

    Args:
        best_predictors_csv:  cfg.output_dir / "chapter2" / "best_predictors.csv"
        experimental_fastas:  list of FASTA paths with experimentally characterised proteins
        output_dir:           cfg.output_dir / "chapter2"
        tmp_dir:              temporary directory for BLAST intermediates
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Building nodes …")
    nodes = build_nodes(best_predictors_csv, experimental_fastas)
    print(f"  {len(nodes)} proteins ({(nodes['source'] == 'PREDICTION').sum()} GWAS, "
          f"{nodes['source'].str.startswith('EXPERIMENTAL').sum()} experimental)")

    print("Running all-vs-all blastp …")
    edges = generate_edges(nodes, tmp_dir)
    non_singleton = edges[edges["query"] != edges["target"]]
    print(f"  {len(non_singleton)} edges after coverage + pident filtering")

    nodes_out = output_dir / "sgnh_network_nodes.tsv"
    edges_out = output_dir / "sgnh_network_edges.tsv"
    nodes.to_csv(nodes_out, sep="\t", index=False)
    edges.to_csv(edges_out, sep="\t", index=False)
    print(f"Saved → {nodes_out.name}, {edges_out.name}")
