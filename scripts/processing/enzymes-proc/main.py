"""
processing/enzymes-proc — experimental enzyme data preparation.

Pipeline:
  1. prepare_prophage_db() → output_dir/manual-outputs/prophage/prophage_db.*
  2. export_sequences()    → enzymes/{proteinID}/sequence.fasta
  3. prepare_af3_json()    → output_dir/manual-upload/enzymes_batch_NNN.json
  4. symlink_structures()  → enzymes/{proteinID}/structure.cif
  5. blast_vs_prophage()   → enzymes/{proteinID}/against-prophages/raw_blast.tsv

Output root: cfg.output_dir / "enzymes"
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from prophage_db import prepare_prophage_db
from enzyme_export import (
    export_sequences,
    prepare_af3_json,
    symlink_structures,
    blast_vs_prophage,
)

cfg = Config()

# ---------------------------------------------------------------------------
# Run flags — toggle steps without modifying config.yml
# When ON: step runs; skips per-protein if output file already exists.
# When OFF: step is skipped entirely.
# ---------------------------------------------------------------------------
RUN_AF3_JSON_PREP  = True   # True → write AF3 JSON files; skips if file exists
RUN_PROPHAGE_BLAST = True   # True → blastp per protein; skips if raw_blast.tsv exists

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
enzymes_xlsx  = cfg.input_dir / "enzymes" / "enzymes.xlsx"
af3_input_dir = cfg.output_dir / "manual-outputs" / "alphafold3" / "enzymes"
enzymes_root  = cfg.output_dir / "enzymes"
blast_db_dir  = cfg.output_dir / "manual-outputs" / "prophage"

# ---------------------------------------------------------------------------
# Step 1 — prophage BLAST DB (checkpoint: skips if DB files already exist)
# ---------------------------------------------------------------------------
print("Step 1: Preparing prophage BLAST DB …")
prepare_prophage_db(
    prophage_faa_pattern = str(cfg.input_dir / "gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa"),
    db_dir               = blast_db_dir,
)

# ---------------------------------------------------------------------------
# Load metadata
# ---------------------------------------------------------------------------
print("Loading enzymes.xlsx …")
enzymes_df = pd.read_excel(enzymes_xlsx, sheet_name="enzymes", engine="openpyxl")
print(f"  {len(enzymes_df)} proteins loaded")

# ---------------------------------------------------------------------------
# Step 2 — sequence.fasta per protein
# ---------------------------------------------------------------------------
print("\nStep 2: Exporting sequences …")
export_sequences(enzymes_df, enzymes_root)

# ---------------------------------------------------------------------------
# Step 3 — AF3 JSON upload files
# ---------------------------------------------------------------------------
if RUN_AF3_JSON_PREP:
    print("\nStep 3: Preparing AF3 JSON upload files …")
    prepare_af3_json(enzymes_df, cfg.output_dir, af3_input_dir)
else:
    print("\n[skip] Step 3: AF3 JSON prep (RUN_AF3_JSON_PREP = False)")

# ---------------------------------------------------------------------------
# Step 4 — symlink AF3 model_0.cif
# ---------------------------------------------------------------------------
print("\nStep 4: Symlinking AF3 structures …")
symlink_structures(enzymes_df, af3_input_dir, enzymes_root)

# ---------------------------------------------------------------------------
# Step 5 — BLASTP vs prophage DB
# ---------------------------------------------------------------------------
print("\nStep 5: BLASTP vs prophage proteins …")
blast_vs_prophage(
    enzymes_df   = enzymes_df,
    enzymes_root = enzymes_root,
    blast_db_dir = blast_db_dir,
    run_blast    = RUN_PROPHAGE_BLAST,
)

print("\nDone.")
