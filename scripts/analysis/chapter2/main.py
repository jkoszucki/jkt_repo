import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from klocus_extraction import KLocusExtractor

cfg = Config()

genbank_dir = cfg.input_dir / "gwas/1_BACTERIA/1_GENBANK_GENOMES"
fasta_dir   = cfg.input_dir / "gwas/1_BACTERIA/2_FASTA_GENOMES_NT"
kloci_dir   = cfg.input_dir / "gwas/4_K_LOCI"
jsonl_path  = cfg.output_dir / "chapter2" / "kaptive_results.jsonl"

extractor = KLocusExtractor(genbank_dir, fasta_dir, kloci_dir, jsonl_path)
extractor.run_kaptive()
extractor.extract()
