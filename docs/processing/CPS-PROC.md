# processing/cps-proc — CPS K-type diversity analysis

Data preparation module feeding **chapter 4**. Parses 81 *Klebsiella pneumoniae* CPS repeating unit structures, extracts modifications, and computes pairwise structural similarity.

## Pipeline

```
cps.xlsx (sheets: ktypes, modifications)
        │
        ├─▶ Step 1: KTypeTablesAPI       → ktypes.csv
        ├─▶ Step 2: KTypeModificationsAPI → ktypes_modifications.csv
        └─▶ Step 3: KTypeSimilarityAPI   → ktypes_sim.csv
```

## Input data

| Path | Description |
|------|-------------|
| `input_dir/cps.xlsx` | CPS structure workbook |

Sheets used:
- `ktypes` — repeating unit structures, bond types, monosaccharide composition (81 K-types: 77 from K-PAM database, 4 from independent literature)
- `modifications` — O-acetyl and O-pyruvyl positions, bond types, and occupancy per K-type

## Output data

Written to `output_dir/chapter3/`:

| File | Description |
|------|-------------|
| `chapter3/ktypes.csv` | Processed K-type table: parsed core/branch linkages, bond types, monosaccharide composition |
| `chapter3/ktypes_modifications.csv` | Modifications per K-type (type, position, bond, occupancy) |
| `chapter3/ktypes_sim.csv` | Pairwise path-based Jaccard similarity (J_core, J_branch, J_total) for all K-type pairs |

## Similarity score definition

Fingerprint = set of all contiguous subpaths of the cyclic core (length 1 to n) plus branch subpaths. Jaccard similarity computed over fingerprint sets. Monosaccharide chirality and chemical modifications deliberately excluded — score captures topology only.

Three scores per pair:
- `J_core` — similarity of core subpaths only
- `J_branch` — similarity of branch subpaths only
- `J_total` — similarity of all subpaths combined

## Library modules (`lib/`)

| Module | Class | Role |
|--------|-------|------|
| `ktypes_base.py` | `BaseKTypeAPI` | XLSX loader base class; sheet caching |
| `ktypes_process.py` | `KTypeTablesAPI` | Parses backbone and branch linkages; writes `ktypes.csv` |
| `ktypes_modifications.py` | `KTypeModificationsAPI` | Exports modifications sheet; writes `ktypes_modifications.csv` |
| `ktypes_similarity.py` | `KTypeSimilarityAPI` | Pairwise fingerprint Jaccard similarity; writes `ktypes_sim.csv` |
| `kpam_get.py` | — | Standalone script to download PDB files from K-PAM (not part of main pipeline) |

## Run

```
cd scripts/processing/cps-proc
conda run -n jkoszucki python main.py
```
