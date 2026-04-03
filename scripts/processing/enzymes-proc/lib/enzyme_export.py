"""
Steps 2–5 of the enzymes-proc pipeline.

  export_sequences()   — sequence.fasta per protein from enzymes.xlsx
  prepare_af3_json()   — AF3 JSON upload files; RBP → trimer, others → monomer
  symlink_structures() — structure.cif symlink → AF3 model_0.cif
  blast_vs_prophage()  — BLASTP each protein vs prophage DB
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _clean_seq(sequence: str) -> str:
    """Strip gaps and whitespace from a sequence string."""
    return sequence.replace("-", "").replace(" ", "").replace("\n", "").replace("\r", "")


def _is_trimer(protein_id: str) -> bool:
    return "RBP" in protein_id


def _make_af3_job(protein_id: str, sequence: str) -> dict:
    count = 3 if _is_trimer(protein_id) else 1
    return {
        "name": protein_id,
        "modelSeeds": [],
        "sequences": [
            {
                "proteinChain": {
                    "sequence": sequence,
                    "count": count,
                }
            }
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }


def _read_seq_index(faa_path: Path) -> dict[str, str]:
    """Return {record_id: sequence} for every record in a FASTA file."""
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(faa_path, "fasta")}


def _run_blastp(
    query_id: str,
    sequence: str,
    db_path: Path,
    evalue: float = 1e-5,
) -> list[tuple]:
    with tempfile.TemporaryDirectory() as tmp:
        q_fasta = Path(tmp) / "query.faa"
        q_fasta.write_text(f">{query_id}\n{sequence}\n")
        result = subprocess.run(
            [
                "blastp",
                "-query", str(q_fasta),
                "-db", str(db_path),
                "-outfmt",
                "6 qseqid sseqid pident length mismatch gapopen "
                "qstart qend sstart send evalue bitscore qcovhsp qlen slen",
                "-evalue", str(evalue),
                "-max_hsps", "1",
            ],
            capture_output=True, text=True, check=True,
        )
    return [
        tuple(line.strip().split("\t"))
        for line in result.stdout.splitlines()
        if line.strip()
    ]


# ---------------------------------------------------------------------------
# Step 2 — sequence.fasta per protein
# ---------------------------------------------------------------------------

def export_sequences(enzymes_df: pd.DataFrame, enzymes_root: Path) -> None:
    """Write enzymes_root/{proteinID_lower}/sequence.fasta for each protein."""
    for _, row in enzymes_df.iterrows():
        protein_id  = str(row["proteinID"])
        protein_dir = enzymes_root / protein_id.lower()
        sequence    = _clean_seq(str(row["sequence"]))
        protein_dir.mkdir(parents=True, exist_ok=True)
        (protein_dir / "sequence.fasta").write_text(f">{protein_id}\n{sequence}\n")
        print(f"  {protein_id.lower()}: sequence.fasta")


# ---------------------------------------------------------------------------
# Step 3 — AF3 JSON upload files
# ---------------------------------------------------------------------------

def prepare_af3_json(
    enzymes_df: pd.DataFrame,
    output_dir: Path,
    af3_input_dir: Path,
    batch_size: int = 30,
) -> None:
    """
    Write one JSON file per batch containing up to batch_size job dicts as a JSON array.
    Skips proteins whose AF3 model_0.cif already exists in af3_input_dir/{proteinID}/.
    Checkpoint: skips a batch file if it already exists.
    With ≤ 30 proteins (minus already downloaded) this produces a single batch_001.json.
    """
    batches_dir = output_dir / "manual-upload"
    batches_dir.mkdir(parents=True, exist_ok=True)

    # Filter out proteins already downloaded
    pending = []
    for _, row in enzymes_df.iterrows():
        protein_id = str(row["proteinID"])
        id_lower   = protein_id.lower()
        cif = af3_input_dir / id_lower / f"fold_{id_lower}_model_0.cif"
        if cif.exists():
            print(f"  [{protein_id}] AF3 structure already downloaded — skipping")
        else:
            pending.append(row)

    if not pending:
        print("  All proteins already downloaded — nothing to write")
        return

    # Split into batches and write
    batches = [pending[i:i + batch_size] for i in range(0, len(pending), batch_size)]

    for batch_idx, batch_rows in enumerate(batches, start=1):
        out = batches_dir / f"enzymes_batch_{batch_idx:03d}.json"

        if out.exists():
            print(f"  [enzymes_batch_{batch_idx:03d}] JSON exists — skipping")
            continue

        jobs = []
        for row in batch_rows:
            protein_id = str(row["proteinID"])
            sequence   = _clean_seq(str(row["sequence"]))
            jobs.append(_make_af3_job(protein_id, sequence))

        out.write_text(json.dumps(jobs, indent=2))
        names = ", ".join(str(r["proteinID"]) for r in batch_rows)
        print(f"  enzymes_batch_{batch_idx:03d}.json — {len(jobs)} proteins: {names}")


# ---------------------------------------------------------------------------
# Step 4 — symlink AF3 model_0.cif
# ---------------------------------------------------------------------------

def symlink_structures(
    enzymes_df: pd.DataFrame,
    af3_input_dir: Path,
    enzymes_root: Path,
) -> None:
    """
    Create enzymes_root/{proteinID_lower}/structure.cif → af3_input_dir/{proteinID_lower}/fold_*_model_0.cif.
    Skips silently if the AF3 folder is not yet present.
    Replaces existing symlink on re-run.
    """
    for _, row in enzymes_df.iterrows():
        protein_id = str(row["proteinID"])
        id_lower   = protein_id.lower()
        src = af3_input_dir / id_lower / f"fold_{id_lower}_model_0.cif"
        if not src.exists():
            print(f"  [{id_lower}] AF3 structure not found — skipping")
            continue

        dest = enzymes_root / id_lower / "structure.cif"
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        dest.symlink_to(src.resolve())
        print(f"  {id_lower}: structure.cif → {src.name}")


# ---------------------------------------------------------------------------
# Step 5 — BLASTP vs prophage DB
# ---------------------------------------------------------------------------

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "qlen", "slen",
]


def blast_vs_prophage(
    enzymes_df: pd.DataFrame,
    enzymes_root: Path,
    blast_db_dir: Path,
    run_blast: bool = True,
) -> None:
    """
    BLASTP each protein against the prophage DB.
    Reads sequence from enzymes_root/{proteinID}/sequence.fasta.
    Writes enzymes_root/{proteinID}/against-prophages/raw_blast.tsv.
    Checkpoint: skips per-protein if raw_blast.tsv already exists.
    """
    db_path  = blast_db_dir / "prophage_db"
    faa_path = blast_db_dir / "prophage_proteins.faa"

    if not list(blast_db_dir.glob("prophage_db.*")):
        print(f"  [error] Prophage BLAST DB not found at {blast_db_dir}")
        print("  Run scripts/helpers/prophage_db.py first.")
        return

    print("  Loading prophage sequence index …")
    seq_index = _read_seq_index(faa_path)
    print(f"  {len(seq_index)} prophage proteins indexed")

    for _, row in enzymes_df.iterrows():
        protein_id = str(row["proteinID"])
        id_lower   = protein_id.lower()
        seq_fasta  = enzymes_root / id_lower / "sequence.fasta"

        if not seq_fasta.exists():
            print(f"  [{id_lower}] sequence.fasta not found — skipping BLAST")
            continue

        against_dir = enzymes_root / id_lower / "against-prophages"
        against_dir.mkdir(parents=True, exist_ok=True)
        raw_tsv = against_dir / "raw_blast.tsv"

        if raw_tsv.exists():
            print(f"  [{id_lower}] raw_blast.tsv exists — skipping")
            continue
        if not run_blast:
            continue

        records = list(SeqIO.parse(seq_fasta, "fasta"))
        if not records:
            print(f"  [{id_lower}] empty sequence.fasta — skipping")
            continue
        sequence = str(records[0].seq)

        print(f"  BLASTing {id_lower} …", flush=True)
        blast_rows = _run_blastp(protein_id, sequence, db_path)
        blast_df   = pd.DataFrame(blast_rows, columns=BLAST_COLS)
        blast_df["sseq"] = blast_df["sseqid"].map(seq_index)
        blast_df.to_csv(raw_tsv, sep="\t", index=False)
        print(f"    {len(blast_df)} hits → raw_blast.tsv")
