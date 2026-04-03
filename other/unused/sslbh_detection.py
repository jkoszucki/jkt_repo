"""
Step 3: Secondary SSLBH (acetyltransferase) detection.

SSLBH (single-stranded left-handed beta-helix, acetyltransferases) cannot be
detected reliably from PC_reported_topology alone. Two methods are used:

Method A — HHsearch (prob ≥ 0.50):
    Re-examine per-PC HHsearch annotations for entries classified as
    other-ecod-reported-topology. If any ECOD hit has SSLBH topology
    with prob ≥ 0.50, classify as SSLBH acetyltransferase.
    Creates symlinks:
        gwas_root/sslbh-ecod-hhsearch/{locus}/{clustering_level}/{PC}/
        → gwas_root/other-ecod-reported-topology/{locus}/{clustering_level}/{PC}/

Method B — TMalign (TMscore ≥ 0.75):
    Compare no-ecod AF3 monomer structures to experimental acetyltransferase
    AF3 structures. Retain entries with TMscore ≥ 0.75 vs any reference.
    Creates symlinks:
        gwas_root/sslbh-tmalign/{locus}/{clustering_level}/{PC}/
        → gwas_root/no-ecod-reported-topology/{locus}/{clustering_level}/{PC}/
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd

SSLBH_PROB_MIN = 0.50
TMALIGN_TM_MIN = 0.75

# Keywords used to identify SSLBH topology in ECOD annotations.
# May need tuning once actual annotation strings are inspected.
_SSLBH_KEYWORDS = ("sslbh", "left-handed beta", "acetyltransfer", "gnat")


def _is_sslbh_topology(topology: str) -> bool:
    """Return True if the topology string indicates an SSLBH acetyltransferase."""
    t = topology.lower()
    return any(kw in t for kw in _SSLBH_KEYWORDS)


def _make_symlink(src: Path, link: Path) -> bool:
    """
    Create a relative symlink link → src.
    Overwrites stale symlinks. Warns and returns False if src does not exist.
    """
    if not src.exists():
        print(f"    [warn] symlink target not found (Step 2 may not have run): {src}")
        return False
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.exists() or link.is_symlink():
        link.unlink()
    rel_target = Path(os.path.relpath(src, link.parent))
    link.symlink_to(rel_target)
    return True


def detect_sslbh_hhsearch(
    gwas_hits_tsv: Path,
    hhsearch_dir: Path,
    gwas_root: Path,
    prob_threshold: float = SSLBH_PROB_MIN,
) -> list[tuple[str, str, str]]:
    """
    Secondary HHsearch lookup for SSLBH among other-ecod entries.

    For each other-ecod (locus, clustering_level, PC) in gwas_hits.tsv:
      - Load per-clustering-level HHsearch annotations
      - Find ECOD hits for this PC with prob >= prob_threshold
      - If any hit has SSLBH topology → create symlink in sslbh-ecod-hhsearch/

    Args:
        gwas_hits_tsv:  gwas_root / "gwas_hits.tsv"
        hhsearch_dir:   input_dir/.../HHSEARCH/ (per-version subdirs)
        gwas_root:      output root; symlinks in sslbh-ecod-hhsearch/
        prob_threshold: minimum HHsearch probability (default 0.50)

    Returns:
        List of (locus, clustering_level, PC) detected as SSLBH.
    """
    df = pd.read_csv(gwas_hits_tsv, sep="\t")
    other = df[df["ecod_type"] == "other-ecod"][
        ["locus", "clustering_level", "PC"]
    ].drop_duplicates()
    print(f"  Checking {len(other)} other-ecod entries for secondary SSLBH …")

    annot_cache: dict[str, pd.DataFrame] = {}
    out_base = gwas_root / "sslbh-ecod-hhsearch"
    detected: list[tuple[str, str, str]] = []

    for _, row in other.iterrows():
        locus, cl, pc = row["locus"], row["clustering_level"], row["PC"]

        if cl not in annot_cache:
            tsv = hhsearch_dir / cl / "clusters_functions.tsv"
            if tsv.exists():
                a = pd.read_csv(tsv, sep="\t")
                if "version" in a.columns:
                    a = a.rename(columns={"version": "clustering_level"})
                annot_cache[cl] = a
            else:
                annot_cache[cl] = pd.DataFrame()

        annot = annot_cache[cl]
        if annot.empty or "PC" not in annot.columns:
            continue

        pc_annot = annot[annot["PC"] == pc]
        if "clustering_level" in pc_annot.columns:
            pc_annot = pc_annot[pc_annot["clustering_level"] == cl]

        if pc_annot.empty or "db" not in pc_annot.columns:
            continue

        ecod_hits = pc_annot[
            (pc_annot["db"] == "ECOD") &
            (pc_annot.get("prob", pd.Series(dtype=float)) >= prob_threshold)
        ]
        if ecod_hits.empty:
            continue

        if "reported_topology_PC" not in ecod_hits.columns:
            continue

        has_sslbh = ecod_hits["reported_topology_PC"].dropna().apply(_is_sslbh_topology).any()
        if not has_sslbh:
            continue

        src  = gwas_root / "other-ecod-reported-topology" / locus / cl / pc
        link = out_base / locus / cl / pc
        if _make_symlink(src, link):
            detected.append((locus, cl, pc))
            print(f"    [sslbh-hhsearch] {locus}/{cl}/{pc}")

    print(f"  {len(detected)} SSLBH detections via HHsearch")
    return detected


def _run_tmalign(query: Path, ref: Path) -> float:
    """Run TMalign; return TM-score normalised by reference length (Chain_2)."""
    result = subprocess.run(
        ["TMalign", str(query), str(ref)],
        capture_output=True, text=True, check=False,
    )
    for line in result.stdout.splitlines():
        if "TM-score=" in line and "Chain_2" in line:
            try:
                return float(line.split("TM-score=")[1].split()[0])
            except (IndexError, ValueError):
                pass
    return 0.0


def detect_sslbh_tmalign(
    gwas_hits_tsv: Path,
    experimental_af3_dir: Path,
    gwas_root: Path,
    tmscore_threshold: float = TMALIGN_TM_MIN,
) -> list[tuple[str, str, str]]:
    """
    TMalign-based SSLBH detection among no-ecod entries.

    For each no-ecod (locus, clustering_level, PC) in gwas_hits.tsv:
      - Locate structure.cif at no-ecod-reported-topology/{locus}/{cl}/{PC}/protein/
      - Run TMalign vs every CIF in experimental_af3_dir/
      - If best TMscore >= tmscore_threshold → create symlink in sslbh-tmalign/

    Args:
        gwas_hits_tsv:          gwas_root / "gwas_hits.tsv"
        experimental_af3_dir:   output_dir/manual-outputs/alphafold3/experimental/
                                (AF3 monomers for experimental acetyltransferases)
        gwas_root:              output root; symlinks in sslbh-tmalign/
        tmscore_threshold:      minimum TMscore (default 0.75)

    Returns:
        List of (locus, clustering_level, PC) detected as SSLBH.
    """
    if shutil.which("TMalign") is None:
        raise RuntimeError(
            "TMalign binary not found. Install TMalign and ensure it is on PATH."
        )

    ref_cifs = sorted(experimental_af3_dir.glob("**/*.cif")) if experimental_af3_dir.exists() else []
    if not ref_cifs:
        print(f"  [warn] No reference CIFs in {experimental_af3_dir} — skipping TMalign")
        return []
    print(f"  {len(ref_cifs)} reference acetyltransferase structures")

    df = pd.read_csv(gwas_hits_tsv, sep="\t")
    no_ecod = df[df["ecod_type"] == "no-ecod"][
        ["locus", "clustering_level", "PC"]
    ].drop_duplicates()
    print(f"  {len(no_ecod)} no-ecod entries for TMalign comparison …")

    out_base = gwas_root / "sslbh-tmalign"
    detected: list[tuple[str, str, str]] = []

    for _, row in no_ecod.iterrows():
        locus, cl, pc = row["locus"], row["clustering_level"], row["PC"]
        query_cif = (
            gwas_root / "no-ecod-reported-topology" / locus / cl / pc
            / "protein" / "structure.cif"
        )
        if not query_cif.exists():
            continue

        best_tm = max((_run_tmalign(query_cif, ref) for ref in ref_cifs), default=0.0)
        if best_tm < tmscore_threshold:
            continue

        src  = gwas_root / "no-ecod-reported-topology" / locus / cl / pc
        link = out_base / locus / cl / pc
        if _make_symlink(src, link):
            detected.append((locus, cl, pc))
            print(f"    [sslbh-tmalign] {locus}/{cl}/{pc}  TM={best_tm:.3f}")

    print(f"  {len(detected)} SSLBH detections via TMalign")
    return detected
