"""
Chapter 2, Figure 1, Panel C
─────────────────────────────
Sequence similarity network for no-ecod representative sequences.

Computes all-vs-all BLASTP of clustered representative sequences from
no-ecod PCs and writes Cytoscape-ready edge and node tables.

BLASTP parameters (per methodology):
    e-value     ≤ 10^-3
    % identity  ≥ 50 %
    query cov   ≥ 80 %
    subject cov ≥ 80 %

Coverage interval labels on edges:
    low     (80–90 %)
    high    (>90 %)

Input:
    sequences_fasta — FASTA of clustered representative sequences for no-ecod
                      PCs.  Placed by the user in output_dir/raw/ after
                      running MMseqs2 easy-cluster on no-ecod representatives.
    node_metadata   — optional CSV with columns [sequence_id, K_locus, PC]
                      to annotate nodes.  If not provided, K_locus and PC are
                      parsed from sequence IDs where possible.

Outputs:
    plots_dir/figure1-panelC_edges.tsv   — Cytoscape edge table
    plots_dir/figure1-panelC_nodes.tsv   — Cytoscape node table
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

# Re-use BLAST machinery from the SGNH network module
from sgnh_network import generate_edges


# ---------------------------------------------------------------------------
# Node assembly
# ---------------------------------------------------------------------------

def _load_nodes(sequences_fasta: Path, node_metadata_csv: Path | None) -> pd.DataFrame:
    """
    Parse a FASTA file into a node DataFrame.
    Columns: query, proteinID, source, seq, K_locus, PC, seq_len
    """
    from Bio import SeqIO

    rows = []
    for rec in SeqIO.parse(sequences_fasta, "fasta"):
        rows.append({
            "proteinID": rec.id,
            "seq":       str(rec.seq),
            "seq_len":   len(rec.seq),
        })
    nodes = pd.DataFrame(rows)

    if node_metadata_csv is not None and node_metadata_csv.exists():
        meta = pd.read_csv(node_metadata_csv)
        nodes = nodes.merge(meta, left_on="proteinID", right_on="sequence_id", how="left")
    else:
        # Try to parse K_locus and PC from typical ID patterns like "KL1_PC0232_repr"
        nodes["K_locus"] = nodes["proteinID"].str.extract(r"(KL\d+)")
        nodes["PC"]      = nodes["proteinID"].str.extract(r"(PC\d+)")

    nodes = nodes.reset_index(drop=True)
    nodes.index = nodes.index + 1
    nodes["source"] = "NO_ECOD_PC"
    nodes["query"]  = nodes.index.astype(str) + "_" + nodes["proteinID"]
    return nodes[["query", "proteinID", "source", "K_locus", "PC", "seq", "seq_len"]]


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def compute_no_ecod_network(
    sequences_fasta: Path,
    plots_dir: Path,
    tmp_dir: Path,
    node_metadata_csv: Path | None = None,
) -> None:
    """
    Build Cytoscape edge and node tables for the no-ecod similarity network.

    Args:
        sequences_fasta:    FASTA with clustered no-ecod representative sequences
        plots_dir:          output directory (scripts/figures/chapter2/plots/)
        tmp_dir:            scratch directory for BLAST intermediates
        node_metadata_csv:  optional CSV with [sequence_id, K_locus, PC] columns
    """
    if not sequences_fasta.exists():
        print(f"  [figure1-panelC] Sequences FASTA not found: {sequences_fasta}\n"
              "  Place the clustered no-ecod representative sequences there and "
              "re-run with PLOT_FIGURE1_PANELC = True")
        return

    print(f"  [figure1-panelC] Loading sequences from {sequences_fasta.name} …")
    nodes = _load_nodes(sequences_fasta, node_metadata_csv)
    print(f"  [figure1-panelC] {len(nodes)} sequences")

    print("  [figure1-panelC] Running all-vs-all BLASTP …")
    edges = generate_edges(
        nodes,
        tmp_dir,
        min_cov=0.80,
        min_pident=50.0,
        max_evalue=1e-3,
        intervals=[80, 90, 100],  # low (80–90%), high (90–100%)
    )

    non_singleton = edges[edges["query"] != edges["target"]]
    print(f"  [figure1-panelC] {len(non_singleton)} edges after filtering")

    # Write Cytoscape tables
    nodes_cols = ["query", "proteinID", "source", "K_locus", "PC", "seq_len"]
    edges_cols = ["query", "target", "pident", "qcov", "scov", "evalue",
                  "bitscore", "coverage_interval", "interaction"]

    nodes_out = plots_dir / "figure1-panelC_nodes.tsv"
    edges_out = plots_dir / "figure1-panelC_edges.tsv"

    nodes[[c for c in nodes_cols if c in nodes.columns]].to_csv(
        nodes_out, sep="\t", index=False
    )
    edges[[c for c in edges_cols if c in edges.columns]].to_csv(
        edges_out, sep="\t", index=False
    )
    print(f"  [figure1-panelC] → {nodes_out.name}, {edges_out.name}")
