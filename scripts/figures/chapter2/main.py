"""
figures/chapter2 — Chapter 2 plots.

Reads from output_dir/gwas-data/ and output_dir/acetyltransferase/.
Writes plots to scripts/figures/chapter2/plots/.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from figure1_panelA import plot_figure1_panelA
from figure2_panelB import plot_figure2_panelB

cfg       = Config()
plots_dir = Path(__file__).resolve().parent / "plots"

gwas_hits_tsv      = cfg.output_dir / "gwas-data" / "gwas_hits.tsv"
pyseer_sslbh_tsv   = cfg.output_dir / "acetyltransferase" / "acetyl-gwas" / "pyseer_hits_sslbh_pvalcor005.tsv"

# ---------------------------------------------------------------------------
# Figure 1A — 6-panel F1 vs MCC scatter coloured by ECOD topology
# ---------------------------------------------------------------------------
print("Figure 1A: GWAS scatter …")
plot_figure1_panelA(
    gwas_hits_tsv = gwas_hits_tsv,
    plots_dir     = plots_dir,
    style         = cfg.style,
)

# ---------------------------------------------------------------------------
# Figure 2B — 6-panel precision vs recall, acetyltransferases highlighted
# ---------------------------------------------------------------------------
print("Figure 2B: precision vs recall scatter …")
plot_figure2_panelB(
    pyseer_tsv = pyseer_sslbh_tsv,
    plots_dir  = plots_dir,
    style      = cfg.style,
)

print("\nDone.")
