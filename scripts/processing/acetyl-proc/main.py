"""
processing/acetyl-proc — acetyltransferase detection and annotation.

Pipeline:
  1. detect_sslbh_pc80()         — detect acetyltransferase PC80s from raw_hhsuite
                                    (ECOD SSLBH / PHROGs keywords / PFAM PF01757)
  2. flag_pyseer_with_sslbh()    — flag GWAS PCs in pyseer_hits_all (pvalue_corr ≤ 0.05)
                                    by tracing PC80 proteins → GWAS clustering levels
                                    → acetyltransferase/acetyl-gwas/pyseer_hits_sslbh_pvalcor005.tsv
  3. export_pc80_annotation()    — write PC80-level annotation tables and per-PC files
                                    → acetyltransferase/acetyl-annot-pc80/
  4. RUN_TMALIGN → subprocess    → acetyl-gwas/no-ecod-reported-topology-tm-align.tsv  (tmtools env)
  5. merge_tmalign_into_pyseer() — pivot TM-align results, merge into pyseer table

Output root: cfg.output_dir / "acetyltransferase"
Feeds figures/chapter2/.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from sslbh_detect import detect_sslbh_pc80
from cluster_map import flag_pyseer_with_sslbh
from pc80_export import export_pc80_annotation
from merge_tmalign import merge_tmalign_into_pyseer
from prophage_db import prepare_prophage_db

cfg = Config()

# ---------------------------------------------------------------------------
# Run flags
# ---------------------------------------------------------------------------
RUN_PROPHAGE_BLAST = True   # True → BLASTP per PC80; skips if raw_blast.tsv exists
RUN_TMALIGN        = True   # True → TM-align no-ecod GWAS predictors vs references (tmtools env)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
raw_hhsuite_tsv        = cfg.input_dir / "gwas/2_PROPHAGES/raw_hhsuite.tsv"
pc2proteins_tsv        = cfg.input_dir / "gwas/2_PROPHAGES/pcs2proteins.tsv"
prophages_metadata_tsv = cfg.input_dir / "gwas/2_PROPHAGES/prophages_metadata.tsv"
pyseer_tsv             = cfg.input_dir / "gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv"
mmseqs_dir             = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/2_MMSEQS"
prophage_faa_pattern   = str(cfg.input_dir / "gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa")

gwas_root        = cfg.output_dir / "gwas-data"
gwas_hits_tsv    = gwas_root / "gwas_hits.tsv"
blast_db_dir     = cfg.output_dir / "manual-outputs" / "prophage"
af3_dir          = cfg.output_dir / "manual-outputs" / "alphafold3"
experimental_dir = af3_dir / "enzymes"

acetyl_dir = cfg.output_dir / "acetyltransferase"

# ---------------------------------------------------------------------------
# Step 1 — Detect SSLBH / acetyltransferase PC80s
# ---------------------------------------------------------------------------
print("Step 1: Detecting SSLBH/acetyltransferase PC80s from raw_hhsuite …")
detected = detect_sslbh_pc80(raw_hhsuite_tsv, pc2proteins_tsv)
print(f"  {len(detected)} PC80s detected")

# Build flat protein → PC80 info map for downstream steps
protein_to_pc80: dict[str, dict] = {}
for _, row in detected.iterrows():
    info = {"pc80": row["pc80"], "detected_by": row["detected_by"]}
    for pid in row["protein_ids"]:
        protein_to_pc80[pid] = info
print(f"  {len(protein_to_pc80)} SSLBH protein IDs indexed")

# ---------------------------------------------------------------------------
# Step 2 — Flag pyseer GWAS hits with SSLBH evidence
# ---------------------------------------------------------------------------
print("\nStep 2: Flagging pyseer GWAS hits with SSLBH evidence …")
flag_pyseer_with_sslbh(
    pyseer_tsv      = pyseer_tsv,
    protein_to_pc80 = protein_to_pc80,
    mmseqs_dir      = mmseqs_dir,
    acetyl_dir      = acetyl_dir,
)

# ---------------------------------------------------------------------------
# Step 3 — Export PC80 annotation
# ---------------------------------------------------------------------------
print("\nStep 3: Exporting PC80 annotation …")
blast_db_path, seq_index = prepare_prophage_db(prophage_faa_pattern, blast_db_dir)
print(f"  {len(seq_index)} prophage proteins in index")

export_pc80_annotation(
    detected               = detected,
    raw_hhsuite_tsv        = raw_hhsuite_tsv,
    pc2proteins_tsv        = pc2proteins_tsv,
    prophages_metadata_tsv = prophages_metadata_tsv,
    seq_index              = seq_index,
    db_path                = blast_db_path,
    af3_dir                = af3_dir,
    acetyl_dir             = acetyl_dir,
    run_blast              = RUN_PROPHAGE_BLAST,
)

# ---------------------------------------------------------------------------
# Step 4 — TM-align screen (tmtools env, via subprocess)
# ---------------------------------------------------------------------------
tmalign_output = acetyl_dir / "acetyl-gwas" / "no-ecod-reported-topology-tm-align.tsv"

if RUN_TMALIGN:
    print("\nStep 4: TM-align screen for no-ecod GWAS predictors …")
    tmalign_script = Path(__file__).resolve().parent / "lib" / "tmalign_screen.py"
    subprocess.run(
        [
            "conda", "run", "-n", "tmtools", "python",
            str(tmalign_script),
            "--gwas-hits",        str(gwas_hits_tsv),
            "--gwas-root",        str(gwas_root),
            "--experimental-dir", str(experimental_dir),
            "--output",           str(tmalign_output),
        ],
        check=True,
    )
else:
    print("\n[skip] Step 4: TM-align screen (RUN_TMALIGN = False)")

# ---------------------------------------------------------------------------
# Step 5 — Merge TM-align results into pyseer table
# ---------------------------------------------------------------------------
print("\nStep 5: Merging TM-align results into pyseer table …")
merge_tmalign_into_pyseer(tmalign_output, acetyl_dir)

print("\nDone.")
