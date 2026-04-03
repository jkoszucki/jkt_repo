"""
Step 7: TM-align screen for no-ecod GWAS predictors vs reference acetyltransferases.

CLI script invoked by main.py via subprocess in the 'tmtools' conda environment:

    conda run -n tmtools python lib/tmalign_screen.py \\
        --gwas-hits  /path/gwas-data/gwas_hits.tsv \\
        --gwas-root  /path/gwas-data/ \\
        --experimental-dir /path/manual-outputs/alphafold3/experimental/ \\
        --output /path/acetyltransferase/tmalign_results.tsv

For each no-ecod GWAS predictor with a structure.cif, compare against all three
reference proteins using tmtools. Both TM-scores (query-normalised and
reference-normalised) and RMSD are recorded.

Output columns:
    locus, clustering_level, pc, reference,
    tm_score_query, tm_score_ref, rmsd
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from Bio.PDB import MMCIFParser
from tmtools.io import get_residue_data
from tmtools import tm_align

REFERENCE_PROTEIN_IDS = [
    "PROTEIN01_MOD_AC_K1",
    "PROTEIN02_MOD_AC_K2",
    "PROTEIN03_MOD_AC_K57",
]

NO_ECOD_FOLDER = "no-ecod-reported-topology"


def run_tmalign_screen(
    gwas_hits_tsv: Path,
    gwas_root: Path,
    experimental_dir: Path,
    output_tsv: Path,
) -> None:
    """
    Screen all no-ecod GWAS predictors against reference acetyltransferases.

    Skips predictors without a structure.cif.
    Appends to output_tsv if it already exists (checkpoint: skips rows already present).
    """
    no_ecod = _read_no_ecod_hits(gwas_hits_tsv)
    print(f"  No-ecod GWAS predictors: {len(no_ecod)}", flush=True)

    # Load reference structures once
    ref_structures = _load_references(experimental_dir)
    if not ref_structures:
        print(
            "  [skip] No reference structures found — expected at "
            f"{experimental_dir}/{{lowercase_protein_id}}/fold_{{lowercase_protein_id}}_model_0.cif",
            flush=True,
        )
        return
    print(f"  Reference structures loaded: {list(ref_structures)}", flush=True)

    # Load existing results to enable checkpointing
    done_keys: set[tuple[str, str, str, str]] = set()
    if output_tsv.exists():
        with open(output_tsv, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                done_keys.add((row["locus"], row["clustering_level"], row["pc"], row["reference"]))
        print(f"  Checkpoint: {len(done_keys)} rows already in {output_tsv.name}", flush=True)

    rows = []
    parser = MMCIFParser(QUIET=True)

    for hit in no_ecod:
        locus = hit["locus"]
        cl    = hit["clustering_level"]
        pc    = hit["PC"]

        cif_path = gwas_root / NO_ECOD_FOLDER / locus / cl / pc / "protein" / "structure.cif"
        if not cif_path.exists():
            print(f"  [skip] {locus}/{cl}/{pc}: structure.cif not found", flush=True)
            continue

        try:
            query_struct = parser.get_structure("query", str(cif_path))
            query_chain  = next(query_struct.get_chains())
            coords_q, seq_q = get_residue_data(query_chain)
        except Exception as exc:
            print(f"  [warn] {locus}/{cl}/{pc}: failed to parse CIF — {exc}", flush=True)
            continue

        for ref_id, (coords_r, seq_r) in ref_structures.items():
            if (locus, cl, pc, ref_id) in done_keys:
                continue
            try:
                result = tm_align(coords_q, coords_r, seq_q, seq_r)
                rows.append({
                    "locus":             locus,
                    "clustering_level":  cl,
                    "pc":                pc,
                    "reference":         ref_id,
                    "tm_score_query":    round(result.tm_norm_chain1, 4),
                    "tm_score_ref":      round(result.tm_norm_chain2, 4),
                    "rmsd":              round(result.rmsd, 4),
                })
                print(
                    f"  {locus}/{cl}/{pc} vs {ref_id}: "
                    f"TM(q)={result.tm_norm_chain1:.3f} "
                    f"TM(r)={result.tm_norm_chain2:.3f} "
                    f"RMSD={result.rmsd:.2f}",
                    flush=True,
                )
            except Exception as exc:
                print(f"  [warn] {locus}/{cl}/{pc} vs {ref_id}: TM-align failed — {exc}", flush=True)

    if not rows:
        print("  No new results to write.", flush=True)
        return

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["locus", "clustering_level", "pc", "reference",
                  "tm_score_query", "tm_score_ref", "rmsd"]
    write_header = not output_tsv.exists()
    with open(output_tsv, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerows(rows)
    total = len(done_keys) + len(rows)
    print(f"  → {output_tsv}  ({total} rows total, {len(rows)} new)", flush=True)


def _read_no_ecod_hits(gwas_hits_tsv: Path) -> list[dict]:
    """Read gwas_hits.tsv and return rows where ecod_type == 'no-ecod'."""
    hits = []
    with open(gwas_hits_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("ecod_type") == "no-ecod":
                hits.append(row)
    return hits


def _load_references(experimental_dir: Path) -> dict[str, tuple]:
    """
    Load CIF structures for all reference proteins.
    Looks for {experimental_dir}/{lowercase_id}/fold_{lowercase_id}_model_0.cif.
    Returns {protein_id: (coords, seq)}.
    """
    parser = MMCIFParser(QUIET=True)
    refs = {}
    for ref_id in REFERENCE_PROTEIN_IDS:
        lid = ref_id.lower()
        cif_path = experimental_dir / lid / f"fold_{lid}_model_0.cif"
        if not cif_path.exists():
            print(f"  [warn] Reference structure not found: {cif_path}", flush=True)
            continue
        try:
            struct = parser.get_structure(ref_id, str(cif_path))
            chain  = next(struct.get_chains())
            coords, seq = get_residue_data(chain)
            refs[ref_id] = (coords, seq)
        except Exception as exc:
            print(f"  [warn] Failed to parse {cif_path}: {exc}", flush=True)
    return refs


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="TM-align screen: no-ecod GWAS vs reference acetyltransferases")
    p.add_argument("--gwas-hits",        required=True, type=Path)
    p.add_argument("--gwas-root",        required=True, type=Path)
    p.add_argument("--experimental-dir", required=True, type=Path)
    p.add_argument("--output",           required=True, type=Path)
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    run_tmalign_screen(
        gwas_hits_tsv    = args.gwas_hits,
        gwas_root        = args.gwas_root,
        experimental_dir = args.experimental_dir,
        output_tsv       = args.output,
    )
