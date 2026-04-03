# processing/enzymes-proc — experimental enzyme data preparation

Data preparation module for experimentally characterised CPS-modifying enzymes. Reads metadata from `enzymes.xlsx`, organises sequences and AF3 structures per protein, and runs each protein against all prophage proteins.

## Pipeline

```
input_dir/gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa
        │
        ▼ Step 1: build prophage BLAST DB (once; checkpoint)
output_dir/manual-outputs/prophage/prophage_db.*

enzymes.xlsx (sheet: enzymes)
        │
        ├─▶ Step 2: export sequence.fasta per protein     → enzymes/{proteinID}/sequence.fasta
        ├─▶ Step 3: prepare AF3 JSON upload files         → output_dir/manual-upload/enzymes_batch_NNN.json (array of ≤30 jobs)
        │
        │   [manual: run AF3, download, place in output_dir/manual-outputs/alphafold3/enzymes/{proteinID_lower}/]
        │
        ├─▶ Step 4: symlink AF3 model_0.cif               → enzymes/{proteinID_lower}/structure.cif → output_dir/manual-outputs/alphafold3/enzymes/{proteinID_lower}/fold_*_model_0.cif
        └─▶ Step 5: BLASTP vs prophage DB                 → enzymes/{proteinID}/against-prophages/raw_blast.tsv
```

## Run flags

| Flag | Default | Controls |
|------|---------|----------|
| `RUN_PROPHAGE_BLAST` | `True` | BLASTP of each protein vs prophage DB; skips per-protein if `raw_blast.tsv` exists |
| `RUN_AF3_JSON_PREP` | `True` | Write AF3 JSON upload files to `output_dir/manual-upload/`; skips if files already exist |

## Input data

| Path | Description |
|------|-------------|
| `input_dir/enzymes/enzymes.xlsx` | Enzyme metadata workbook (sheet: `enzymes`) |
| `output_dir/manual-outputs/alphafold3/enzymes/{proteinID_lower}/` | AF3 raw output folder per protein, placed manually after download |
| `input_dir/gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa` | Prophage protein FASTAs (shared with `gwas-proc`) |

### `enzymes.xlsx` — sheet `enzymes`

Columns used:

| Column | Description |
|--------|-------------|
| `proteinID` | Unique protein identifier — used as folder name and FASTA header |
| `sequence` | Amino acid sequence |
| `K-type` | Associated K-type (may be empty for unknown) |
| `source` | Organism source (`KP`, `AB`, `EC`) |
| `modification` | Enzymatic activity (`O-acetylation`, `O-pyruvylation`, `deacetylation`) |

### Proteins (10 total)

| proteinID | K-type | source | modification | AF3 run-as |
|-----------|--------|--------|--------------|------------|
| PROTEIN01_MOD_AC_K1 | K1 | KP | O-acetylation | monomer |
| PROTEIN02_MOD_AC_K2 | K2 | KP | O-acetylation | monomer |
| PROTEIN03_MOD_AC_K57 | K57 | KP | O-acetylation | monomer |
| PROTEIN04_MOD_PYR_K1 | K1 | KP | O-pyruvylation | monomer |
| PROTEIN05_MOD_DAC_K2 | K2 | KP | deacetylation | monomer |
| PROTEIN06_MOD_DAC_UNK | — | KP | deacetylation | monomer |
| PROTEIN07_RBP_DAC_K2 | K2 | KP | deacetylation | **trimer** |
| PROTEIN08_RBP_DAC_K26 | K26 | AB | deacetylation | **trimer** |
| PROTEIN09_RBP_DAC_K55 | K55 | KP | deacetylation | **trimer** |
| PROTEIN10_RBP_DAC_LPS4s | LPS4s | EC | deacetylation | **trimer** |

**Run-as rule:** proteins with `RBP` in `proteinID` → homotrimer; all others → monomer.

## AlphaFold3 conventions

JSON upload files are prepared by `prepare_af3_json()` in `lib/enzyme_export.py`. Use `proteinID` as the AF3 job name. Each batch file is a JSON array of up to 30 job dicts, written to `output_dir/manual-upload/enzymes_batch_NNN.json`. With 10 proteins this produces a single `enzymes_batch_001.json`.

After downloading from the AF3 server, place each raw output folder at:
```
output_dir/manual-outputs/alphafold3/enzymes/{proteinID_lower}/
    fold_{proteinid_lower}_model_0.cif   ← used by Step 4
    fold_{proteinid_lower}_model_1.cif
    ...
    fold_{proteinid_lower}_full_data_0.json
    ...
```

Step 4 creates a symlink `output_dir/enzymes/{proteinID_lower}/structure.cif` → `output_dir/manual-outputs/alphafold3/enzymes/{proteinID_lower}/fold_{proteinid_lower}_model_0.cif`. Skips silently for any protein whose AF3 folder has not yet been placed.

## Output data

Written to `output_dir/enzymes/`:

```
enzymes/
    {proteinID_lower}/
        sequence.fasta              — single-record FASTA; header = proteinID; sequence from enzymes.xlsx
        structure.cif               — symlink → output_dir/manual-outputs/alphafold3/enzymes/{proteinID_lower}/fold_*_model_0.cif
        against-prophages/
            raw_blast.tsv           — BLASTP hits vs all prophage proteins
                                      columns: qseqid, sseqid, pident, length, mismatch, gapopen,
                                               qstart, qend, sstart, send, evalue, bitscore,
                                               qcovhsp, qlen, slen, sseq
```

## Prophage BLAST DB

Built by `prepare_prophage_db()` from `scripts/helpers/prophage_db.py`, called as Step 1 of `main.py`. Shared with `gwas-proc`. Checkpoint: if `prophage_db.*` already exist, Step 1 returns immediately.

What it does:
1. Reads all `input_dir/gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa`; extracts protein name from `' ||| '`-delimited headers (`>CDS_1 ||| KPN_B1_PHAGE002_M_PROTEIN_1 ||| ...`)
2. Writes `prophage_proteins.faa` (all proteins, clean `>{protein_name}` headers)
3. Builds `prophage_db.*` from `prophage_proteins.faa`

```
output_dir/manual-outputs/prophage/
    prophage_proteins.faa   — all prophage proteins, curated headers
    prophage_db.*           — blastp DB files (makeblastdb)
```

For BLAST column definitions see the output structure section above.

## Library modules (`lib/`)

| Module | Function | Role |
|--------|----------|------|
| `helpers/prophage_db.py` | `prepare_prophage_db()` | Step 1: build prophage FASTA + BLAST DB; checkpoint on existing DB files |
| `enzyme_export.py` | `export_sequences()` | Step 2: write `sequence.fasta` per protein from `enzymes.xlsx` |
| `enzyme_export.py` | `prepare_af3_json()` | Step 3: write AF3 JSON upload files; RBP → trimer, others → monomer; checkpoint per file |
| `enzyme_export.py` | `symlink_structures()` | Step 4: create `structure.cif` symlink → AF3 `model_0.cif`; skips missing |
| `enzyme_export.py` | `blast_vs_prophage()` | Step 5: BLASTP each protein vs prophage DB; writes `against-prophages/raw_blast.tsv`; checkpoint on file existence |

## Run

```
cd scripts/processing/enzymes-proc
conda run -n jkoszucki python main.py
```
