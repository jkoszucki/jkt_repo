"""
Prepare AlphaFold3 JSON upload files from a list of sequence FASTA paths.

Usage
-----
Pass a flat list of paths to sequence.fasta files — GWAS predictors,
experimental proteins, or any mix.  The function determines IDs and metadata
from the path structure alone; no TSV input is required.

GWAS paths
----------
Paths matching the GWAS output hierarchy:

    {ecod_type}/{locus}/{clustering_level}/{PC}/protein/sequence.fasta

are assigned an internal ID that encodes (clustering_level, PC):

    PCI50C50 + PC001  →  pci50c50_pc001
    PCI80C50 + PC001  →  pci80c50_pc001   ← different cluster despite same label

This solves the duplicate-PC-name problem: PC001 at PCI50C50 and PC001 at
PCI80C50 are different clusters and get different IDs.

When the same (clustering_level, PC) is the best predictor for multiple
K-loci, only one AF3 job is submitted; af3_placement.tsv records every
locus that should receive a copy of the structure.

Non-GWAS paths
--------------
Paths that do not match the GWAS pattern (experimental proteins, etc.) are
assigned sequential IDs:  PROTEIN001, PROTEIN002, ...

Outputs written to output_dir/manual-upload/
--------------------------------------------
gwas_batch_001.json       — JSON array of ≤30 AF3 job dicts
gwas_batch_002.json
...
af3_manifest.tsv          — af3_id | run_as | batch | sequence
af3_placement.tsv         — af3_id | source_path | locus | ecod_type
                            | clustering_level | PC | protein_id

The placement table is read by the organise/symlink step after structures
are downloaded from the AF3 server.

AF3 conventions (from CLAUDE.md)
---------------------------------
- model_0 is always used when organising downloaded results.
- sgnh-ecod-reported-topology and ssrbh-ecod-reported-topology → trimer
  (chains A, B, C, same sequence).
- All other ecod types → monomer (chain A only).
- Batch size ≤ 30  (AF3 server daily submission limit).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# ECOD types predicted as trimers (SGNH hydrolases and SSRBH depolymerases).
_TRIMER_ECOD_TYPES = {
    "sgnh-ecod-reported-topology",
    "ssrbh-ecod-reported-topology",
}

# Matches the GWAS output hierarchy inside a path string.
# Anchored on the known clustering_level format (PCI<id>C<cov>) and PC<n>.
# Groups: ecod_type, locus, clustering_level, PC.
_GWAS_PATH_RE = re.compile(
    r"(?P<ecod_type>[^/]+)"
    r"/(?P<locus>[^/]+)"
    r"/(?P<clustering_level>PCI\d+C\d+)"
    r"/(?P<PC>PC\d+)"
    r"/protein/sequence\.fasta$"
)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_gwas_path(path: Path) -> dict | None:
    """Extract GWAS metadata from path, or None if the path does not match."""
    m = _GWAS_PATH_RE.search(path.as_posix())
    return m.groupdict() if m else None


def _make_af3_id(clustering_level: str, pc: str) -> str:
    """Unique internal ID encoding (clustering_level, PC).

    Examples:
        PCI50C50, PC001  →  pci50c50_pc001
        PCI80C50, PC001  →  pci80c50_pc001
    """
    return f"{clustering_level.lower()}_{pc.lower()}"


def _read_first_sequence(path: Path) -> tuple[str, str] | None:
    """Return (protein_id, sequence) for the first record in a FASTA file.

    Gaps are stripped so the sequence is suitable for AF3 input.
    Returns None when the file is missing or empty.
    """
    if not path.exists():
        return None
    for record in SeqIO.parse(path, "fasta"):
        return record.id, str(record.seq).replace("-", "")
    return None


def _make_job_dict(af3_id: str, sequence: str, run_as: str) -> dict:
    """Return the AF3 JSON job dict for one prediction.

    Trimers use count=3; monomers use count=1.
    """
    count = 3 if run_as == "trimer" else 1
    return {
        "name": af3_id,
        "modelSeeds": [],
        "sequences": [
            {
                "proteinChain": {
                    "sequence": sequence,
                    "count": count,
                }
            }
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def prepare_af3_batches(
    sequence_paths: list[Path],
    output_dir: Path,
    batch_size: int = 30,
) -> None:
    """Prepare AF3 JSON upload files and companion manifest tables.

    Args:
        sequence_paths: paths to sequence.fasta files — GWAS predictors and/or
                        experimental proteins, in any combination and order.
        output_dir:     destination for manual-upload/ folder and manifest files.
        batch_size:     maximum predictions per batch folder (default 30,
                        matching the AF3 server daily submission limit).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    batches_dir = output_dir / "manual-upload"

    # ------------------------------------------------------------------
    # Pass 1: parse every path, read sequences, assign IDs.
    # ------------------------------------------------------------------
    # GWAS deduplication: (clustering_level, PC) → af3_id already assigned.
    gwas_seen: dict[tuple[str, str], str] = {}
    # Sequential counter for non-GWAS proteins.
    protein_counter = 0

    placement_rows: list[dict] = []
    # Ordered dict preserves insertion order for stable batch assignment.
    # af3_id → (sequence, run_as, ecod_type)
    predictions: dict[str, tuple[str, str, str]] = {}
    missing: list[str] = []

    for path in sequence_paths:
        result = _read_first_sequence(path)
        if result is None:
            missing.append(str(path))
            continue
        protein_id, sequence = result

        meta = _parse_gwas_path(path)

        if meta is not None:
            # ---- GWAS path: ID encodes (clustering_level, PC) ----
            cl  = meta["clustering_level"]
            pc  = meta["PC"]
            key = (cl, pc)

            if key not in gwas_seen:
                af3_id = _make_af3_id(cl, pc)
                run_as = "trimer" if meta["ecod_type"] in _TRIMER_ECOD_TYPES else "monomer"
                gwas_seen[key] = af3_id
                predictions[af3_id] = (sequence, run_as, meta["ecod_type"])
            else:
                af3_id = gwas_seen[key]

            placement_rows.append({
                "af3_id":           af3_id,
                "source_path":      str(path),
                "locus":            meta["locus"],
                "ecod_type":        meta["ecod_type"],
                "clustering_level": cl,
                "PC":               pc,
                "protein_id":       protein_id,
            })

        else:
            # ---- Non-GWAS path: sequential PROTEIN id ----
            protein_counter += 1
            af3_id = f"PROTEIN{protein_counter:03d}"
            predictions[af3_id] = (sequence, "monomer", "")

            placement_rows.append({
                "af3_id":           af3_id,
                "source_path":      str(path),
                "locus":            "",
                "ecod_type":        "",
                "clustering_level": "",
                "PC":               "",
                "protein_id":       protein_id,
            })

    if missing:
        print(f"  [warn] {len(missing)} paths missing or empty — skipped:")
        for m in missing:
            print(f"    {m}")

    # ------------------------------------------------------------------
    # Pass 2: group by ECOD type; batch each group independently.
    # Files named {ecod_type}-gwas-batch_{NNN}.json — no mixing across types.
    # ------------------------------------------------------------------
    batches_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict] = []

    # Group predictions by ecod_type (stable order within each group)
    by_ecod: dict[str, list[tuple[str, str, str]]] = {}
    for af3_id, (sequence, run_as, ecod_type) in predictions.items():
        by_ecod.setdefault(ecod_type, []).append((af3_id, sequence, run_as))

    total_batches = 0
    for ecod_type in sorted(by_ecod):
        prefix = f"{ecod_type}-gwas" if ecod_type else "gwas"
        group  = by_ecod[ecod_type]
        chunks = [group[i:i + batch_size] for i in range(0, len(group), batch_size)]
        for batch_num, chunk in enumerate(chunks, start=1):
            filename = f"{prefix}-batch_{batch_num:03d}.json"
            out  = batches_dir / filename
            jobs = []
            for af3_id, sequence, run_as in chunk:
                jobs.append(_make_job_dict(af3_id, sequence, run_as))
                manifest_rows.append({
                    "af3_id":    af3_id,
                    "run_as":    run_as,
                    "ecod_type": ecod_type,
                    "batch":     filename,
                    "sequence":  sequence,
                })
            out.write_text(json.dumps(jobs, indent=2))
            print(f"  {filename} — {len(jobs)} predictions")
            total_batches += 1

    # ------------------------------------------------------------------
    # Write tables.
    # ------------------------------------------------------------------
    pd.DataFrame(placement_rows).to_csv(batches_dir / "af3_placement.tsv", sep="\t", index=False)
    pd.DataFrame(manifest_rows).to_csv(batches_dir / "af3_manifest.tsv", sep="\t", index=False)

    n_pred   = len(predictions)
    n_placed = len(placement_rows)
    print(f"  {n_placed} input paths → {n_pred} unique predictions")
    print(f"  {n_pred} predictions → {total_batches} batch file(s) of ≤{batch_size} across {len(by_ecod)} ECOD type(s)")
    print(f"  Batches written to: {batches_dir}/")
