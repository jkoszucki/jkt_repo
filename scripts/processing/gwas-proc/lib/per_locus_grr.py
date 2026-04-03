"""
Compute wGRR similarity between the reference K-locus and all category
genomes (match, locus_only, pc_only) for each best-predictor row.

Uses GRRpair methodology (Kupczok et al. 2022 / de Sousa et al. 2021):
  BBH via BLAST (evalue < 0.1) + powerneedle + wGRR = Σ(ident/100)/min(g1,g2)

Output per predictor:
    gwas_root/{ecod_folder}/{locus}/{clustering_level}/{PC}/locus/grr.tsv
        genomeID  category  wgrr  grr_sum  genes_ref  genes_genome  n_bbh
"""

from __future__ import annotations

import re
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO

EVALUE     = 0.1
CATEGORIES = ["match", "locus_only", "pc_only"]
REF_ID     = "REFERENCE"


# ---------------------------------------------------------------------------
# Protein extraction
# ---------------------------------------------------------------------------

def _sanitize(s: str) -> str:
    """Replace characters that break BLAST/powerneedle sequence IDs."""
    return re.sub(r"[\[\]|.\s]", "_", s)


def _extract_proteins(gb_path: Path, sample_id: str) -> list[tuple[str, str]]:
    """Return [(label, seq), ...] where label = sample_id__protid."""
    records = []
    try:
        for rec in SeqIO.parse(gb_path, "genbank"):
            rec_id = _sanitize(rec.id)
            for feat in rec.features:
                if feat.type != "CDS" or "translation" not in feat.qualifiers:
                    continue
                prot_id = f"{rec_id}_{int(feat.location.start)}_{int(feat.location.end)}"
                label   = f"{sample_id}__{prot_id}"
                records.append((label, feat.qualifiers["translation"][0]))
    except Exception as exc:
        print(f"    [warn] could not parse {gb_path.name}: {exc}")
    return records


def _write_faa(proteins: list[tuple[str, str]], path: Path) -> None:
    with open(path, "w") as fh:
        for label, seq in proteins:
            fh.write(f">{label}\n{seq}\n")


# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------

def _run_blast(merged_faa: Path, db: Path, out: Path) -> None:
    subprocess.run(
        ["blastp", "-query", str(merged_faa), "-db", str(db),
         "-outfmt", "7", "-evalue", str(EVALUE), "-out", str(out)],
        check=True, capture_output=True,
    )


# ---------------------------------------------------------------------------
# BBH extraction (same logic as klocus_grr.py, filtered to reference pairs)
# ---------------------------------------------------------------------------

def _parse_bbh(
    blast_tsv: Path,
    genome_ids: set[str],
) -> dict[str, list[tuple[str, str]]]:
    """
    Extract BBH pairs between REFERENCE and each genome in genome_ids.

    Returns {genome_id: [(ref_bare_prot, genome_bare_prot), ...]}
    The bare protein IDs are what comes after '__' and are used directly
    in the powerneedle pair FAA and .bbh files.
    """
    best: dict[tuple[str, str, str], str] = {}
    with open(blast_tsv) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 11:
                continue
            qid, sid = parts[0], parts[1]
            try:
                evalue = float(parts[10])
            except ValueError:
                continue
            if evalue >= EVALUE:
                continue
            q_sample, q_prot = qid.split("__", 1)
            s_sample, s_prot = sid.split("__", 1)
            if q_sample == s_sample:
                continue
            # Only track pairs that involve REFERENCE
            if q_sample != REF_ID and s_sample != REF_ID:
                continue
            key = (q_sample, q_prot, s_sample)
            if key not in best:
                best[key] = s_prot

    bbh: dict[str, list[tuple[str, str]]] = {}
    seen: set[tuple[str, str, str, str]] = set()
    for (q_s, q_p, s_s), s_p in best.items():
        if best.get((s_s, s_p, q_s)) != q_p:
            continue
        ref_s, gen_s = (q_s, s_s) if q_s == REF_ID else (s_s, q_s)
        ref_p, gen_p = (q_p, s_p) if q_s == REF_ID else (s_p, q_p)
        if gen_s not in genome_ids:
            continue
        key4 = (ref_p, gen_s, gen_p)
        if key4 not in seen:
            seen.add(key4)
            bbh.setdefault(gen_s, []).append((ref_p, gen_p))

    return bbh


# ---------------------------------------------------------------------------
# powerneedle + wGRR
# ---------------------------------------------------------------------------

def _compute_wgrr(
    genome_id: str,
    bbh_pairs: list[tuple[str, str]],
    all_proteins: dict[str, str],
    tmp_dir: Path,
) -> tuple[float, float, int, int, int]:
    """
    Run powerneedle for one reference–genome pair.
    Returns (wgrr, grr_sum, genes_ref, genes_genome, n_bbh).
    """
    ref_bare    = {k.split("__", 1)[1]: v for k, v in all_proteins.items()
                   if k.startswith(REF_ID + "__")}
    genome_bare = {k.split("__", 1)[1]: v for k, v in all_proteins.items()
                   if k.startswith(genome_id + "__")}

    g_ref    = len(ref_bare)
    g_genome = len(genome_bare)

    base  = tmp_dir / genome_id[:60]   # truncate for filesystem safety
    faa_p = base.with_suffix(".faa")
    bbh_p = base.with_suffix(".bbh")
    ndl_p = base.with_suffix(".needle")
    alg_p = base.with_suffix(".needlealg")

    # Pair FAA uses bare protein IDs (no sample prefix) — as in GRRpair
    _write_faa(list(ref_bare.items()) + list(genome_bare.items()), faa_p)
    bbh_p.write_text("".join(f"{r_p}\t{g_p}\n" for r_p, g_p in bbh_pairs))

    subprocess.run(
        ["powerneedle", "-gapopen", "10", "-gapextend", "0.5", "-brief",
         "-pairs", str(bbh_p), "-identities", str(ndl_p),
         "-alignment", str(alg_p), str(faa_p)],
        capture_output=True,
    )

    grr_sum = 0.0
    if ndl_p.exists():
        for line in ndl_p.read_text().splitlines():
            if not line.startswith("Sequence1"):
                try:
                    grr_sum += float(line.split()[2]) / 100
                except (IndexError, ValueError):
                    pass

    for p in (faa_p, bbh_p, ndl_p, alg_p):
        p.unlink(missing_ok=True)

    denom = min(g_ref, g_genome)
    wgrr  = grr_sum / denom if denom > 0 else 0.0
    return wgrr, grr_sum, g_ref, g_genome, len(bbh_pairs)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def _pc_dir(gwas_root: Path, row) -> Path:
    """Return the 4D base directory for a gwas_hits row."""
    return gwas_root / row["ecod_folder"] / row["locus"] / row["clustering_level"] / row["PC"]


def export_grr_vs_reference(
    gwas_hits_tsv: Path,
    gwas_root: Path,
) -> None:
    """
    For each row in gwas_hits.tsv, compute wGRR between the reference GenBank
    and every genome in match/, locus_only/, pc_only/.

    Output: gwas_root/{ecod_folder}/{locus}/{clustering_level}/{PC}/locus/grr.tsv
        genomeID  category  wgrr  grr_sum  genes_ref  genes_genome  n_bbh
    """
    best = pd.read_csv(gwas_hits_tsv, sep="\t")

    for _, row in best.iterrows():
        label     = f"{row['ecod_folder']}/{row['locus']}/{row['clustering_level']}/{row['PC']}"
        locus_dir = _pc_dir(gwas_root, row) / "locus"
        ref_gb    = locus_dir / "reference.gb"
        out_tsv   = locus_dir / "grr.tsv"

        if out_tsv.exists():
            print(f"  [{label}] grr.tsv exists, skipping")
            continue

        if not ref_gb.exists():
            print(f"  [{label}] reference GB not found, skipping")
            continue

        # Collect category genomes
        category_genomes: list[tuple[str, str, Path]] = []  # (genomeID, category, gb_path)
        for cat in CATEGORIES:
            cat_dir = locus_dir / cat
            if not cat_dir.exists():
                continue
            for gb in sorted(cat_dir.glob("*.gb")):
                category_genomes.append((gb.stem, cat, gb))

        if not category_genomes:
            print(f"  [{label}] no category genomes found, skipping")
            continue

        n_total = len(category_genomes)
        print(f"  [{label}] {n_total} genomes …", flush=True)

        # Extract proteins: REFERENCE + all category genomes
        ref_prots = _extract_proteins(ref_gb, REF_ID)
        if not ref_prots:
            print(f"  [{label}] reference has no CDS, skipping")
            continue

        all_proteins: dict[str, str] = {label: seq for label, seq in ref_prots}
        genome_ids: set[str] = set()
        for genome_id, _cat, gb_path in category_genomes:
            for label, seq in _extract_proteins(gb_path, genome_id):
                all_proteins[label] = seq
            genome_ids.add(genome_id)

        # BLAST + BBH + powerneedle in a temp directory
        with tempfile.TemporaryDirectory() as _tmp:
            tmp_dir    = Path(_tmp)
            merged_faa = tmp_dir / "merged.faa"
            _write_faa(list(all_proteins.items()), merged_faa)

            db_path = tmp_dir / "blastdb" / "locus_db"
            db_path.parent.mkdir()
            subprocess.run(
                ["makeblastdb", "-in", str(merged_faa), "-dbtype", "prot",
                 "-out", str(db_path)],
                check=True, capture_output=True,
            )

            blast_out = tmp_dir / "blast.tsv"
            _run_blast(merged_faa, db_path, blast_out)

            bbh      = _parse_bbh(blast_out, genome_ids)
            ndl_tmp  = tmp_dir / "needle"
            ndl_tmp.mkdir()

            rows = []
            for genome_id, cat, _ in category_genomes:
                pairs = bbh.get(genome_id, [])
                if not pairs:
                    g_genome = sum(
                        1 for k in all_proteins if k.startswith(genome_id + "__")
                    )
                    rows.append({
                        "genomeID":     genome_id, "category":     cat,
                        "wgrr":         0.0,       "grr_sum":      0.0,
                        "genes_ref":    len(ref_prots), "genes_genome": g_genome,
                        "n_bbh":        0,
                    })
                    continue
                wgrr, grr_sum, g_ref, g_genome, n_bbh = _compute_wgrr(
                    genome_id, pairs, all_proteins, ndl_tmp
                )
                rows.append({
                    "genomeID":     genome_id,          "category":     cat,
                    "wgrr":         round(wgrr, 6),     "grr_sum":      round(grr_sum, 4),
                    "genes_ref":    g_ref,              "genes_genome": g_genome,
                    "n_bbh":        n_bbh,
                })

        result_df = pd.DataFrame(rows, columns=[
            "genomeID", "category", "wgrr", "grr_sum", "genes_ref", "genes_genome", "n_bbh"
        ])
        result_df.to_csv(out_tsv, sep="\t", index=False)

        counts = result_df.groupby("category").size().to_dict()
        print(f"    " + ", ".join(f"{c}={counts.get(c,0)}" for c in CATEGORIES)
              + f" → {out_tsv.name}")
