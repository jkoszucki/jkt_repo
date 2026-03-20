# Chapter 2 — Analysis

Processes the CPS (capsule polysaccharide) structure dataset and exports three CSV tables used downstream by `scripts/figures/chapter2/`.

## Input

| File | Path |
|---|---|
| CPS structure workbook | `input_dir/cps.xlsx` |

Sheets used: `ktypes`, `modifications`.

## Outputs

Written to `output_dir/analysis/chapter2/`:

| File | Description |
|---|---|
| `ktypes.csv` | Processed K-type table with parsed core/branch fields, bond types, and monosaccharide composition |
| `ktypes_modifications.csv` | Raw modifications table (OAc, OPy, etc.) per K-type |
| `ktypes_sim.csv` | Pairwise Jaccard similarity scores across all K-type pairs |

## Run

```
conda run -n jkoszucki python scripts/analysis/chapter2/main.py
```

## Library modules (`lib/`)

| Module | Class | Role |
|---|---|---|
| `ktypes_base.py` | `BaseKTypeAPI` | XLSX loader base class; sheet caching |
| `ktypes_process.py` | `KTypeTablesAPI` | Parses backbone and branch linkages; computes Jaccard similarity |
| `ktypes_modifications.py` | `KTypeModificationsAPI` | Exports the modifications sheet |
| `ktypes_similarity.py` | `KTypeSimilarityAPI` | Thin wrapper for similarity export |
| `kpam_get.py` | — | Standalone script to download PDB files from K-PAM |
