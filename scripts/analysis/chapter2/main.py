import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from gwas_predictor_selection import select_best_predictors
from klocus_grr import compute_klocus_grr

cfg = Config()

table_path  = cfg.input_dir / "gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv"
mmseqs_dir  = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/2_MMSEQS"
k_loci_dir  = cfg.input_dir / "gwas/4_K_LOCI"
output_dir  = cfg.output_dir / "chapter2"

select_best_predictors(table_path, mmseqs_dir, output_dir)
compute_klocus_grr(k_loci_dir, output_dir)
