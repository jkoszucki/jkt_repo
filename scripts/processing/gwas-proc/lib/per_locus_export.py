"""
Step 2: Export per-PC data into the 4D output hierarchy.

For each hit in gwas_hits.tsv:
    gwas_root/{ecod_folder}/{locus}/{clustering_level}/{PC}/
        protein/
            pc.fasta               — MMseqs2 alignment FASTA (full cluster)
            sequence.fasta         — representative sequence (first from alignment)
            structure.cif          — AF3 model_0 (trimer for sgnh/ssrbh, monomer otherwise)
            against-prophages/
                raw_blast.tsv      — BLASTP hits vs all prophage proteins; columns:
                                     qseqid, sseqid, pident, length, mismatch, gapopen,
                                     qstart, qend, sstart, send, evalue, bitscore,
                                     qcovhsp, qlen, slen, sseq

        locus/
            reference.gb           — reference K-locus GenBank
            match/                 — K-locus and PC co-occur
            locus_only/            — K-locus present, PC absent
            pc_only/               — PC present, K-locus not typed


"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

from prophage_db import prepare_prophage_db

import pandas as pd


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def _pc_dir(gwas_root: Path, row) -> Path:
    """Return the 4D base directory for a gwas_hits row."""
    return gwas_root / row["ecod_folder"] / row["locus"] / row["clustering_level"] / row["PC"]


def _read_first_fasta(fasta_path: Path) -> tuple[str, str]:
    """Return (protein_id, sequence) for the first record in a FASTA file."""
    protein_id, seq_lines, found = "", [], False
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if found:
                    break
                protein_id = line[1:]
                found = True
            elif found:
                seq_lines.append(line)
    return protein_id, "".join(seq_lines)



# ---------------------------------------------------------------------------
# AF3 indexing
# ---------------------------------------------------------------------------

def _index_af3(af3_dir: Path) -> dict[str, Path]:
    """
    Walk af3_dir recursively and return {lowercase_key: prediction_folder_path}.

    A prediction folder is any directory that directly contains a
    fold_*_model_0.cif file. The search descends up to 4 levels to handle
    arbitrary batch/grouping subfolder nesting (e.g.
    af3_dir/no-ecod-reported-topology/folds_2026_.../pci80c80_pc0001/).

    Canonical name wins over duplicate runs (e.g. pc001 over pc001_2).
    """
    available: dict[str, Path] = {}
    if not af3_dir.exists():
        print(f"  [warn] AF3 dir not found: {af3_dir} — skipping structure copy")
        return available

    def _add(d: Path) -> None:
        key = d.name.lower()
        if re.search(r"_\d+$", key):
            canonical = re.sub(r"_\d+$", "", key)
            if canonical in available:
                print(f"  [skip duplicate] {d.name}")
                return
        available[key] = d

    # Find all prediction folders (contain fold_*_model_0.cif) up to depth 4
    for cif in sorted(af3_dir.glob("**/fold_*_model_0.cif")):
        _add(cif.parent)

    return available


def _lookup_af3(available: dict[str, Path], *keys: str) -> Path | None:
    """Try each key in order; return the first match or None."""
    for key in keys:
        if key in available:
            return available[key]
    return None


def _copy_cif(folder_path: Path, dest: Path) -> bool:
    """Copy model_0 CIF from an AF3 prediction folder to dest. Returns True on success."""
    cif_src = folder_path / f"fold_{folder_path.name}_model_0.cif"
    if not cif_src.exists():
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(cif_src, dest)
    return True


# ---------------------------------------------------------------------------
# Public: alignment FASTAs and AF3 structures
# ---------------------------------------------------------------------------

def export_per_locus(
    gwas_hits_tsv: Path,
    mmseqs_dir: Path,
    af3_dir: Path,
    gwas_root: Path,
) -> None:
    """
    Copy alignment FASTAs (→ pc.fasta) and AF3 CIF structures (→ structure.cif)
    into the 4D output hierarchy.

    AF3 lookup order per predictor:
      1. {clustering_level}_{PC} (disambiguated af3_id for duplicate PC names)
      2. {PC} (simple name)

    Args:
        gwas_hits_tsv: gwas_root / "gwas_hits.tsv"
        mmseqs_dir:    input_dir/.../2_MMSEQS
        af3_dir:       output_dir/manual-outputs/alphafold3/
        gwas_root:     output root
    """
    hits = pd.read_csv(gwas_hits_tsv, sep="\t")

    # --- Alignments → pc.fasta ---
    print("  Alignments (pc.fasta) …")
    align_ok, align_missing = 0, []
    for _, row in hits.iterrows():
        src = mmseqs_dir / row["clustering_level"] / "alignments" / f"{row['PC']}.fasta"
        if not src.exists():
            align_missing.append(f"{row['locus']}/{row['clustering_level']}/{row['PC']}")
            continue
        dest_dir = _pc_dir(gwas_root, row) / "protein"
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest_dir / "pc.fasta")
        align_ok += 1
    print(f"    Copied {align_ok}/{len(hits)}")
    for m in align_missing[:10]:
        print(f"    missing: {m}")
    if len(align_missing) > 10:
        print(f"    … and {len(align_missing) - 10} more")

    # --- Index AF3 folders ---
    available = _index_af3(af3_dir)

    # --- CIFs → structure.cif ---
    print("  AF3 structures (structure.cif) …")
    cif_ok, cif_missing = 0, []

    for _, row in hits.iterrows():
        pc, cl = row["PC"], row["clustering_level"]
        folder_path = _lookup_af3(
            available,
            f"{cl}_{pc}".lower(),
            pc.lower(),
        )
        dest_dir = _pc_dir(gwas_root, row) / "protein"
        if folder_path and _copy_cif(folder_path, dest_dir / "structure.cif"):
            cif_ok += 1
        else:
            rel = f"{row['ecod_folder']}/{row['locus']}/{cl}/{pc}"
            cif_missing.append((f"{cl}_{pc}", rel))

    total = len(hits)
    print(f"    Copied {cif_ok}/{total}")
    if cif_missing:
        print(f"    Missing ({len(cif_missing)}) — resubmit to AF3:")
        for job_name, rel in cif_missing[:20]:
            print(f"      {job_name:<35} → {rel}/structure.cif")
        if len(cif_missing) > 20:
            print(f"      … and {len(cif_missing) - 20} more")


# ---------------------------------------------------------------------------
# Public: sequence FASTAs
# ---------------------------------------------------------------------------

def export_per_protein(
    gwas_hits_tsv: Path,
    mmseqs_dir: Path,
    gwas_root: Path,
) -> None:
    """
    Write sequence.fasta for each predictor.
    Representative sequence = first record from the MMseqs2 alignment FASTA.
    → gwas_root/{ecod_folder}/{locus}/{cl}/{PC}/protein/sequence.fasta
    """
    hits = pd.read_csv(gwas_hits_tsv, sep="\t")

    for _, row in hits.iterrows():
        locus, pc, cl = row["locus"], row["PC"], row["clustering_level"]
        fasta_path = mmseqs_dir / cl / "alignments" / f"{pc}.fasta"
        if not fasta_path.exists():
            print(f"  [warn] {locus}/{cl}/{pc}: alignment FASTA not found")
            continue
        protein_id, seq = _read_first_fasta(fasta_path)
        if not seq:
            print(f"  [warn] {locus}/{cl}/{pc}: empty sequence")
            continue

        protein_dir = _pc_dir(gwas_root, row) / "protein"
        protein_dir.mkdir(parents=True, exist_ok=True)
        header = protein_id if protein_id else f"{locus}_{pc}"
        with open(protein_dir / "sequence.fasta", "w") as fh:
            fh.write(f">{header}\n{seq}\n")


# ---------------------------------------------------------------------------
# Public: K-locus GenBank dataset
# ---------------------------------------------------------------------------

KPSC_SPECIES = {
    "Klebsiella pneumoniae",
    "Klebsiella variicola subsp. variicola",
    "Klebsiella quasipneumoniae subsp. quasipneumoniae",
    "Klebsiella quasipneumoniae subsp. similipneumoniae",
}

GOOD_CONF = {"Perfect", "Very high", "High", "Good"}


def export_klocus_genbank_dataset(
    gwas_hits_tsv: Path,
    bacteria_metadata_tsv: Path,
    k_loci_dir: Path,
    ref_dir: Path,
    mmseqs_dir: Path,
    gwas_root: Path,
) -> None:
    """
    For each hit in gwas_hits.tsv, build the locus/ dataset in the 4D hierarchy:
        gwas_root/{ecod_folder}/{locus}/{clustering_level}/{PC}/locus/
            reference.gb
            match/        — K-locus and PC co-occur
            locus_only/   — K-locus present, PC absent
            pc_only/      — PC present, K-locus not typed
    """
    hits = pd.read_csv(gwas_hits_tsv, sep="\t")

    meta_raw = pd.read_csv(bacteria_metadata_tsv, sep="\t")
    meta = meta_raw[meta_raw["species"].isin(KPSC_SPECIES)][
        ["genomeID", "K_locus", "K_type", "K_locus_confidence"]
    ].copy()
    print(f"  Loaded metadata: {len(meta)} Kpsc genomes")

    # Fallback index for non-standard filenames (e.g. AB371295.gb for K62)
    # Maps K_type stem (e.g. "K62") or K_locus (e.g. "KL62") → file path
    ref_fallback: dict[str, Path] = {}
    for gb in sorted(ref_dir.glob("*.gb")):
        if re.match(r"^KL?\d+\.gb$", gb.name):
            continue
        try:
            content = gb.read_text(errors="ignore")
            m = re.search(r'K locus[:\s]+KL(\d+)', content)
            if m:
                ref_fallback[f"KL{m.group(1)}"] = gb
                continue
            m = re.search(r'serotype="K(\d+)"', content)
            if m:
                ref_fallback[f"K{m.group(1)}"] = gb
        except Exception:
            pass
    if ref_fallback:
        print(f"  Fallback ref index: {len(ref_fallback)} non-standard files mapped")

    binary_cache: dict[str, pd.DataFrame] = {}

    for _, row in hits.iterrows():
        locus, pc, cl = row["locus"], row["PC"], row["clustering_level"]
        locus_out = _pc_dir(gwas_root, row) / "locus"
        locus_out.mkdir(parents=True, exist_ok=True)

        # Reference GenBank
        locus_meta = meta[meta["K_locus"] == locus]
        k_type_series = locus_meta["K_type"].dropna()
        if not k_type_series.empty:
            k_type = k_type_series.iloc[0]
            m = re.match(r"unknown \((KL\d+)\)", k_type)
            ref_stem = m.group(1) if m else k_type
            ref_src = ref_dir / f"{ref_stem}.gb"
            if not ref_src.exists():
                ref_src = ref_fallback.get(ref_stem) or ref_fallback.get(locus)
            if ref_src and ref_src.exists():
                shutil.copy2(ref_src, locus_out / "reference.gb")
            else:
                print(f"  [warn] {locus}: reference not found for {ref_stem}")
        else:
            print(f"  [warn] {locus}: no K_type in metadata")

        # Binary matrix (cached by clustering level)
        if cl not in binary_cache:
            bm_path = mmseqs_dir / cl / "3_binary_matrix.tsv"
            if bm_path.exists():
                binary_cache[cl] = pd.read_csv(bm_path, sep="\t", index_col=0)
            else:
                print(f"  [warn] binary matrix not found: {bm_path}")
                binary_cache[cl] = pd.DataFrame()

        bm = binary_cache[cl]
        if bm.empty or pc not in bm.index:
            pc_presence = pd.Series(dtype=float)
        else:
            pc_presence = bm.loc[pc]

        if locus_meta.empty:
            print(f"  [warn] {locus}: no Kpsc genomes in metadata")
            continue

        all_df = meta.copy()
        all_df["pc_present"] = all_df["genomeID"].map(pc_presence).fillna(0).astype(int)

        good_locus_mask = (
            (all_df["K_locus"] == locus) &
            all_df["K_locus_confidence"].isin(GOOD_CONF)
        )
        match_df      = all_df[good_locus_mask & (all_df["pc_present"] == 1)]
        locus_only_df = all_df[good_locus_mask & (all_df["pc_present"] == 0)]
        pc_only_df    = all_df[~good_locus_mask & (all_df["pc_present"] == 1)]

        for group_name, group_df in [
            ("match",      match_df),
            ("locus_only", locus_only_df),
            ("pc_only",    pc_only_df),
        ]:
            group_dir = locus_out / group_name
            group_dir.mkdir(exist_ok=True)
            for genome_id in group_df["genomeID"]:
                src = k_loci_dir / f"{genome_id}.gb"
                if src.exists():
                    shutil.copy2(src, group_dir / f"{genome_id}.gb")

        print(
            f"  {locus}/{row['ecod_folder']}/{cl}/{pc}: "
            f"match={len(match_df)}, locus_only={len(locus_only_df)}, pc_only={len(pc_only_df)}"
        )


# ---------------------------------------------------------------------------
# Public: BLAST vs prophage proteins
# ---------------------------------------------------------------------------


def _run_blastp(
    query_id: str, query_seq: str, db_path: Path, tmp_dir: Path, evalue: float = 1e-5
) -> list[tuple[str, ...]]:
    """Run blastp; return list of raw hit rows (outfmt 6 fields)."""
    q_fasta = tmp_dir / f"{query_id}_q.faa"
    q_fasta.write_text(f">{query_id}\n{query_seq}\n")
    result = subprocess.run(
        ["blastp", "-query", str(q_fasta), "-db", str(db_path),
         "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen",
         "-evalue", str(evalue), "-max_hsps", "1"],
        capture_output=True, text=True, check=True,
    )
    return [
        tuple(line.strip().split("\t"))
        for line in result.stdout.splitlines()
        if line.strip()
    ]


def _blast_one_protein(
    query_id: str,
    seq: str,
    db_path: Path,
    seq_index: dict[str, str],
    db_dir: Path,
    protein_dir: Path,
    run_blast: bool = True,
) -> None:
    """
    BLAST one protein against the prophage DB.
    Writes protein_dir/against-prophages/raw_blast.tsv with subject sequences appended.
    Skips if raw_blast.tsv already exists.
    """
    blast_cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "qlen", "slen",
    ]
    against_dir = protein_dir / "against-prophages"
    against_dir.mkdir(parents=True, exist_ok=True)

    raw_tsv = against_dir / "raw_blast.tsv"
    if raw_tsv.exists():
        print(f"  [{query_id}] raw_blast.tsv exists, skipping", flush=True)
        return
    if not run_blast:
        return

    print(f"  BLASTing {query_id} …", flush=True)
    blast_rows = _run_blastp(query_id, seq, db_path, db_dir)
    blast_df = pd.DataFrame(blast_rows, columns=blast_cols)
    blast_df["sseq"] = blast_df["sseqid"].map(seq_index)
    blast_df.to_csv(raw_tsv, sep="\t", index=False)
    print(f"    {len(blast_df)} hits → {raw_tsv.name}", flush=True)


def blast_repr_vs_prophage(
    gwas_hits_tsv: Path,
    mmseqs_dir: Path,
    prophage_faa_pattern: str,
    gwas_root: Path,
    blast_db_dir: Path,
    run_blast: bool = True,
    max_proteins: int | None = None,
) -> None:
    """
    BLAST all predictor proteins against prophage proteins.
    Representative sequence per predictor = first record from alignment FASTA.

    Args:
        max_proteins: if set, stop after this many proteins are BLASTed (for testing).

    Output: gwas_root/{ecod_folder}/{locus}/{cl}/{PC}/protein/against-prophages/raw_blast.tsv
    """
    print("  Prophage BLAST DB …")
    db_path, seq_index = prepare_prophage_db(prophage_faa_pattern, blast_db_dir)
    print(f"  {len(seq_index)} prophage proteins indexed")

    n_blasted = 0
    hits = pd.read_csv(gwas_hits_tsv, sep="\t")
    for _, row in hits.iterrows():
        if max_proteins is not None and n_blasted >= max_proteins:
            break
        locus, pc, cl = row["locus"], row["PC"], row["clustering_level"]
        fasta_path = mmseqs_dir / cl / "alignments" / f"{pc}.fasta"
        if not fasta_path.exists():
            print(f"  [{locus}/{cl}/{pc}] alignment FASTA not found — skipping BLAST")
            continue
        protein_id, rep_seq = _read_first_fasta(fasta_path)
        if not rep_seq:
            print(f"  [{locus}/{cl}/{pc}] empty sequence — skipping")
            continue
        _blast_one_protein(
            query_id=protein_id or f"{locus}_{pc}", seq=rep_seq,
            db_path=db_path, seq_index=seq_index, db_dir=blast_db_dir,
            protein_dir=_pc_dir(gwas_root, row) / "protein",
            run_blast=run_blast,
        )
        n_blasted += 1
