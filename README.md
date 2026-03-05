## Overview
This repository contains code used to process the data and generate all figures used in the PhD thesis titled "XXX" by Janusz Koszucki. Code is organised to reflect scripts that are used in each of the corresponding chapters.

## Repository Structure

```
.
├── config
│   └── config.yml          # Master configuration file
├── README.md
└── scripts
    ├── analysis            # Data processing scripts
    │   └── chapter2
    │       ├── main.py     # Entry point: loads Config, runs analysis pipeline
    │       └── lib/        # Importable modules
    │           ├── ktypes_base.py          # BaseKTypeAPI — XLSX loader base class
    │           ├── kpam_get.py             # KPAM PDB downloader
    │           ├── ktypes_process.py       # KTypeTablesAPI — main processing
    │           ├── ktypes_modifications.py # KTypeModificationsAPI
    │           └── ktypes_similarity.py    # KTypeSimilarityAPI
    ├── figures             # Scripts generating figures
    │   └── chapter2
    │       ├── main.py     # Entry point: loads Config, produces all figures
    │       ├── lib/        # Importable modules (no analysis/ imports)
    │       │   ├── ktypes_draw.py          # KTypeStructureDrawer
    │       │   ├── ktypes_plots.py         # KTypePlotAPI — reads CSVs
    │       │   ├── branch_core_heatmap.py  # Hexbin similarity heatmap
    │       │   ├── composition_groups.py   # Monosaccharide composition groups
    │       │   └── mw_hist.py              # Molecular weight histogram
    │       ├── plots/      # Generated plots (written by main.py)
    │       ├── figures/    # Manual workspace (not touched by scripts)
    │       ├── cps_draw/   # Static structure assets
    │       ├── ktypes_network.cys
    │       └── modifying_enzymes.pdf
    └── helpers
        └── config.py       # Config class — loads config/config.yml
```

## Config

`Config` is a Python class in `scripts/helpers/config.py` that loads `config/config.yml` and exposes paths and style settings as attributes. It is imported once in each `main.py` and values are passed as arguments to all library modules.

```python
from config import Config
cfg = Config()
# cfg.input_dir  — raw input data root
# cfg.output_dir — analysis CSV output root (chapter subdirs live directly here)
# cfg.style      — Style dataclass with font sizes, DPI, and colours
```

Master configuration file (`config/config.yml`):

| Key | Description |
|---|---|
| `paths.input_dir` | Root for raw input data (e.g. `cps.xlsx`) |
| `paths.output_dir` | Root for analysis CSV outputs; chapter subdirs created here |
| `style.*` | Shared plot style — fonts, DPI, colours (see Style section below) |

## Style

Shared visual style is defined once in `config/config.yml` under the `style:` key and loaded by `Config` into a `Style` dataclass (`cfg.style`). All figure scripts receive it as `style=cfg.style`.

```yaml
style:
  dpi: 200
  axis_label_fontsize: 12
  axis_label_fontweight: bold
  tick_fontsize: 8
  tick_fontweight: bold
  kpam_color: "#1f77b4"
  lit_color: "#d62728"
  pyruvylation_color: "#4f81bd"
  acetylation_color: "#f2c94c"
  gray_color: "#bfbfbf"
  both_color: "#DC143C"
```

The default values follow the style of `year.png` (`plot_cumulative_structures`). Add or update keys in `config.yml` to change the style globally across all figures.

## Environment

All scripts run in the `jkoszucki` conda environment (Python 3.10):

```
conda run -n jkoszucki python scripts/analysis/chapter2/main.py
```

## Running chapter 2

**Step 1 — Analysis** (produces 3 CSVs in `data/output/analysis/chapter2/`):
```
conda run -n jkoszucki python scripts/analysis/chapter2/main.py
```
Reads: `data/input/cps.xlsx`
Outputs: `ktypes.csv`, `ktypes_modifications.csv`, `ktypes_sim.csv`

**Step 2 — Figures** (reads CSVs, writes plots to `scripts/figures/chapter2/plots/`):
```
conda run -n jkoszucki python scripts/figures/chapter2/main.py
```

The two steps are independent: figures can be re-run any time after analysis has produced its CSVs.
