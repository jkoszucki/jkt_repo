"""
Export per-locus data files: alignment FASTAs and AF3 ranked_0 CIF structures.

For each best predictor:
    output_dir/{locus}/protein/{locus}_{PC}.fasta   — MMseqs2 alignment
    output_dir/{locus}/protein/{locus}_{PC}.cif     — AF3 model_0 structure (when available)

For each experimental protein (see _EXPERIMENTAL_INFO):
    output_dir/experimental/{subfolder}/{stem}.cif

Multiple AF3 batch folders inside raw_af3_dir are supported.
Duplicate runs (e.g. pc1246_2) are skipped in favour of the canonical folder.
"""

from __future__ import annotations

import glob as _glob
import io
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def _load_experimental_ids(fasta_paths: list[Path]) -> list[str]:
    """Return protein IDs (first token of each FASTA header) from all FASTAs."""
    ids = []
    for fasta in fasta_paths:
        if not fasta.exists():
            continue
        with open(fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    ids.append(line[1:].split()[0].strip())
    return ids


def _experimental_placement(job_name: str) -> tuple[str, str]:
    """
    Return (folder, stem) for an experimental protein AF3 job name.

    164_08_KL2       → (KL2,          KL2_164_08)
    174_38_KL55      → (KL55,         KL55_174_38)
    QNO11465DOT1_K26 → (experimental, K26_QNO11465DOT1)   # NCBI accession → experimental
    YP_004782195DOT1 → (experimental, YP_004782195DOT1)
    """
    m = re.search(r"_(KL?\d+)$", job_name, re.IGNORECASE)
    if m:
        locus   = m.group(1).upper()
        protein = job_name[: m.start()]
        stem    = f"{locus}_{protein}"
        if "DOT" in protein:          # NCBI accession number → experimental
            return "experimental", stem
        return locus, stem
    return "experimental", job_name


def _index_raw_af3(raw_af3_dir: Path) -> dict[str, Path]:
    """
    Walk all batch subdirs in raw_af3_dir and return
    {lowercase_protein_key: protein_folder_path}.
    Canonical name wins over duplicate runs (e.g. pc1246 over pc1246_2).
    """
    available: dict[str, Path] = {}
    batch_dirs = sorted(d for d in raw_af3_dir.iterdir() if d.is_dir())
    if not batch_dirs:
        print(f"  [warn] No batch folders found in {raw_af3_dir}")
        return available

    print(f"  AF3 batch folders ({len(batch_dirs)}): "
          + ", ".join(b.name for b in batch_dirs))

    for batch in batch_dirs:
        for d in batch.iterdir():
            if not d.is_dir():
                continue
            key = d.name.lower()
            if re.search(r"_\d+$", key):
                canonical = re.sub(r"_\d+$", "", key)
                if canonical in available:
                    print(f"  [skip duplicate] {d.name}")
                    continue
            available[key] = d

    return available


# ---------------------------------------------------------------------------
# Experimental protein registry
# ---------------------------------------------------------------------------

# subfolder  : directory name under per-locus/experimental/
# stem       : prefix for all output files in that folder
# cif_folder : where export_per_locus copied the AF3 CIF (relative to per_locus_dir)
# cif_stem   : AF3 CIF filename stem (before .cif)
_EXPERIMENTAL_INFO: dict[str, dict[str, str]] = {
    "164_08_KL2": {
        "subfolder":  "164_08_KL2_4NPA",
        "stem":       "164_08_KL2",
        "cif_folder": "KL2",
        "cif_stem":   "KL2_164_08",
    },
    "174_38_KL55": {
        "subfolder":  "174_38_KL55_4NPA",
        "stem":       "174_38_KL55",
        "cif_folder": "KL55",
        "cif_stem":   "KL55_174_38",
    },
    "QNO11465.1_K26": {
        "subfolder":  "QNO11465.1_K26",
        "stem":       "QNO11465.1_K26",
        "cif_folder": "experimental",
        "cif_stem":   "K26_QNO11465DOT1",
    },
    "YP_004782195.1": {
        "subfolder":  "YP_004782195.1_OAg_LPS",
        "stem":       "YP_004782195.1",
        "cif_folder": "experimental",
        "cif_stem":   "YP_004782195DOT1",
    },
}


def export_per_locus(
    best_predictors_csv: Path,
    mmseqs_dir: Path,
    raw_af3_dir: Path,
    experimental_fastas: list[Path],
    output_dir: Path,
) -> None:
    """
    Copy alignment FASTAs and AF3 CIF structures into per-locus subfolders.

    Args:
        best_predictors_csv:  cfg.output_dir / "chapter2" / "best_predictors.csv"
        mmseqs_dir:           cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/2_MMSEQS"
        raw_af3_dir:          cfg.output_dir / "chapter2" / "raw_af3"
        experimental_fastas:  list of paths to experimental protein FASTAs
        output_dir:           cfg.output_dir / "chapter2" / "per-locus"
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    best = pd.read_csv(best_predictors_csv)

    # --- Step 1: alignments ---
    print("  Alignments ...")
    align_ok, align_missing = 0, []
    for _, row in best.iterrows():
        src = mmseqs_dir / row["version"] / "alignments" / f"{row['PC']}.fasta"
        if not src.exists():
            align_missing.append(f"{row['locus']} ({src.name})")
            continue
        dest_dir = output_dir / row["locus"] / "protein"
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest_dir / f"{row['locus']}_{row['PC']}.fasta")
        align_ok += 1
    print(f"    Copied {align_ok}/{len(best)}")
    for m in align_missing:
        print(f"    missing: {m}")

    # --- Step 2: index AF3 folders ---
    available = _index_raw_af3(raw_af3_dir)

    # --- Step 3: build job map ---
    # Predictor entries: key=PC, value=(locus, stem)
    # Experimental entries: key=job_name, value=(protein_id,) — resolved via _EXPERIMENTAL_INFO
    predictor_map: dict[str, tuple[str, str]] = {}  # PC → (locus, stem)
    for _, row in best.iterrows():
        predictor_map[row["PC"]] = (row["locus"], f"{row['locus']}_{row['PC']}")

    exp_id_map: dict[str, str] = {}  # job_name → protein_id
    for protein_id in _load_experimental_ids(experimental_fastas):
        job_name = protein_id.replace(".", "DOT")
        exp_id_map[job_name] = protein_id

    # --- Step 4: copy CIFs ---
    print("  AF3 structures ...")
    cif_ok, cif_missing = 0, []

    # Predictors → {locus}/protein/
    for pc, (locus, stem) in sorted(predictor_map.items()):
        folder_path = available.get(pc.lower())
        if folder_path is None:
            cif_missing.append((pc, f"{locus}/protein", stem))
            continue
        cif_src = folder_path / f"fold_{folder_path.name}_model_0.cif"
        if not cif_src.exists():
            cif_missing.append((pc, f"{locus}/protein", stem))
            continue
        dest_dir = output_dir / locus / "protein"
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(cif_src, dest_dir / f"{stem}.cif")
        cif_ok += 1

    # Experimental → experimental/{subfolder}/
    for job_name, protein_id in sorted(exp_id_map.items()):
        info = _EXPERIMENTAL_INFO.get(protein_id)
        if info is None:
            continue
        folder_path = available.get(job_name.lower())
        if folder_path is None:
            cif_missing.append((job_name, f"experimental/{info['subfolder']}", info["stem"]))
            continue
        cif_src = folder_path / f"fold_{folder_path.name}_model_0.cif"
        if not cif_src.exists():
            cif_missing.append((job_name, f"experimental/{info['subfolder']}", info["stem"]))
            continue
        dest_dir = output_dir / "experimental" / info["subfolder"]
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(cif_src, dest_dir / f"{info['stem']}.cif")
        cif_ok += 1

    print(f"    Copied {cif_ok}/{len(predictor_map) + len(exp_id_map)}")
    if cif_missing:
        print(f"    Missing ({len(cif_missing)}) — resubmit to AF3:")
        for job_name, folder, stem in cif_missing:
            print(f"      {job_name:<30} → {folder}/{stem}.cif")


KPSC_SPECIES = {
    "Klebsiella pneumoniae",
    "Klebsiella variicola subsp. variicola",
    "Klebsiella quasipneumoniae subsp. quasipneumoniae",
    "Klebsiella quasipneumoniae subsp. similipneumoniae",
}

GOOD_CONF = {"Perfect", "Very high", "High", "Good"}


def export_klocus_genbank_dataset(
    best_predictors_csv: Path,
    bacteria_metadata_tsv: Path,
    k_loci_dir: Path,
    ref_dir: Path,
    mmseqs_dir: Path,
    output_dir: Path,
) -> None:
    """
    For each best-predictor K-locus, build a structured dataset under
    output_dir/{KL}/ with reference GenBank, per-genome GenBanks split by
    PC presence and confidence, and a GRR similarity fragment.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Load shared data ---
    best = pd.read_csv(best_predictors_csv)

    meta_raw = pd.read_csv(bacteria_metadata_tsv, sep="\t")
    meta = meta_raw[meta_raw["species"].isin(KPSC_SPECIES)][
        ["genomeID", "K_locus", "K_type", "K_locus_confidence"]
    ].copy()
    print(f"  Loaded metadata: {len(meta)} Kpsc genomes")

    binary_cache: dict[str, pd.DataFrame] = {}

    # --- Per-locus loop ---
    for _, row in best.iterrows():
        locus = row["locus"]
        pc = row["PC"]
        version = row["version"]

        locus_out = output_dir / locus / "locus"
        locus_out.mkdir(parents=True, exist_ok=True)

        # a. Reference GenBank
        locus_meta = meta[meta["K_locus"] == locus]
        k_type_series = locus_meta["K_type"].dropna()
        if k_type_series.empty:
            print(f"  [warn] {locus}: no K_type found in metadata, skipping reference")
        else:
            k_type = k_type_series.iloc[0]
            m = re.match(r"unknown \((KL\d+)\)", k_type)
            ref_stem = m.group(1) if m else k_type
            ref_src = ref_dir / f"{ref_stem}.gb"
            ref_dest = locus_out / f"{locus}_reference.gb"
            if ref_src.exists():
                shutil.copy2(ref_src, ref_dest)
            else:
                print(f"  [warn] {locus}: reference not found at {ref_src}")

        # b. Load binary matrix (cached)
        if version not in binary_cache:
            bm_path = mmseqs_dir / version / "3_binary_matrix.tsv"
            if bm_path.exists():
                binary_cache[version] = pd.read_csv(bm_path, sep="\t", index_col=0)
            else:
                print(f"  [warn] binary matrix not found: {bm_path}")
                binary_cache[version] = pd.DataFrame()

        bm = binary_cache[version]
        if bm.empty or pc not in bm.index:
            pc_presence = pd.Series(dtype=float)
        else:
            pc_presence = bm.loc[pc]

        # c. Classify Kpsc genomes
        if locus_meta.empty:
            print(f"  [warn] {locus}: no Kpsc genomes found, skipping")
            continue

        all_df = meta.copy()
        all_df["pc_present"] = all_df["genomeID"].map(pc_presence).fillna(0).astype(int)

        good_locus_mask = (all_df["K_locus"] == locus) & all_df["K_locus_confidence"].isin(GOOD_CONF)
        match_df      = all_df[good_locus_mask & (all_df["pc_present"] == 1)]
        locus_only_df = all_df[good_locus_mask & (all_df["pc_present"] == 0)]
        pc_only_df    = all_df[~good_locus_mask & (all_df["pc_present"] == 1)]

        # d. Copy GenBank files
        groups = [
            ("match",      match_df),
            ("locus_only", locus_only_df),
            ("pc_only",    pc_only_df),
        ]
        for group_name, group_df in groups:
            group_dir = locus_out / group_name
            group_dir.mkdir(exist_ok=True)
            for genome_id in group_df["genomeID"]:
                src = k_loci_dir / f"{genome_id}.gb"
                if src.exists():
                    shutil.copy2(src, group_dir / f"{genome_id}.gb")
                else:
                    print(f"  [warn] {locus}/{group_name}: missing {src.name}")

        print(f"  {locus}: match={len(match_df)}, locus_only={len(locus_only_df)}, pc_only={len(pc_only_df)}")


# ---------------------------------------------------------------------------
# Prophage BLAST helpers
# ---------------------------------------------------------------------------

def _build_blast_db(prophage_faa_pattern: str, db_dir: Path) -> tuple[Path, dict[str, str]]:
    """Concatenate prophage FASTAs, build blastp DB. Returns (db_path, seq_index).

    Each sequence is labelled {prophage_stem}__{protein_id} so BLAST sseqid
    encodes both the source prophage and the protein within it.
    """
    seq_index: dict[str, str] = {}
    combined = db_dir / "prophage_combined.faa"
    with open(combined, "w") as out:
        for faa in sorted(_glob.glob(prophage_faa_pattern)):
            prophage = Path(faa).stem
            for rec in SeqIO.parse(faa, "fasta"):
                label = f"{prophage}__{rec.id}"
                seq_index[label] = str(rec.seq)
                out.write(f">{label}\n{rec.seq}\n")
    if not seq_index:
        raise FileNotFoundError(f"No sequences matched: {prophage_faa_pattern}")
    db_path = db_dir / "prophage_db"
    subprocess.run(
        ["makeblastdb", "-in", str(combined), "-dbtype", "prot", "-out", str(db_path)],
        check=True, capture_output=True,
    )
    return db_path, seq_index


def _run_blastp(query_id: str, query_seq: str, db_path: Path, tmp_dir: Path,
                evalue: float = 1e-5) -> list[tuple[str, ...]]:
    """Run blastp, return list of raw hit rows (outfmt 6 fields)."""
    q_fasta = tmp_dir / f"{query_id}_q.faa"
    q_fasta.write_text(f">{query_id}\n{query_seq}\n")
    result = subprocess.run(
        ["blastp", "-query", str(q_fasta), "-db", str(db_path),
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore qcovs",
         "-evalue", str(evalue), "-max_hsps", "1"],
        capture_output=True, text=True, check=True,
    )
    rows = []
    for line in result.stdout.splitlines():
        if line.strip():
            rows.append(tuple(line.strip().split("\t")))
    return rows


def _write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    """Write (label, seq) pairs to a FASTA file."""
    with open(path, "w") as fh:
        for label, seq in records:
            fh.write(f">{label}\n{seq}\n")


def _align_fasta(fasta_path: Path, out_path: Path) -> None:
    """Run MAFFT on fasta_path, write aligned output to out_path."""
    result = subprocess.run(
        ["mafft", "--auto", "--quiet", str(fasta_path)],
        capture_output=True, text=True, check=True,
    )
    out_path.write_text(result.stdout)


def _load_experimental_records(
    fastas: list[Path],
) -> list[tuple[str, dict[str, str], str]]:
    """
    Load experimental protein sequences from FASTA files.
    Returns list of (protein_id, info_dict, seq) where info_dict comes from
    _EXPERIMENTAL_INFO and contains subfolder, stem, cif_folder, cif_stem.
    Proteins not listed in _EXPERIMENTAL_INFO are skipped with a warning.
    """
    records = []
    for fasta in fastas:
        if not fasta.exists():
            continue
        for rec in SeqIO.parse(fasta, "fasta"):
            info = _EXPERIMENTAL_INFO.get(rec.id)
            if info is None:
                print(f"  [warn] experimental protein {rec.id!r} not in _EXPERIMENTAL_INFO — skipping")
                continue
            records.append((rec.id, info, str(rec.seq)))
    return records


def _blast_one_protein(
    query_id: str,
    stem: str,
    seq: str,
    db_path: Path,
    seq_index: dict[str, str],
    db_dir: Path,
    protein_dir: Path,
    bitscore_threshold: int,
    run_blast: bool = True,
) -> None:
    """
    BLAST one protein against the prophage DB and write all outputs to protein_dir:
        {stem}_blast_raw.tsv
        {stem}_blast_filtered.tsv
        {stem}_unaligned.fasta
        {stem}_prophage.fasta

    run_blast=True  — always run blastp and overwrite raw TSV, then filter + align.
    run_blast=False — skip blastp; use existing raw TSV if present, then filter + align.
                      If no raw TSV exists, warn and skip entirely.
    """
    blast_cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qcovs"]
    protein_dir.mkdir(parents=True, exist_ok=True)

    raw_tsv = protein_dir / f"{stem}_blast_raw.tsv"
    if run_blast:
        print(f"  BLASTing {query_id} …", flush=True)
        blast_rows = _run_blastp(query_id, seq, db_path, db_dir)
        blast_df = pd.DataFrame(blast_rows, columns=blast_cols)
        blast_df["bitscore"] = pd.to_numeric(blast_df["bitscore"], errors="coerce")
        blast_df.to_csv(raw_tsv, sep="\t", index=False)
    else:
        if not raw_tsv.exists():
            print(f"  [{query_id}] no raw BLAST found and run_blast=False — skipping", flush=True)
            return
        print(f"  [{query_id}] using existing raw BLAST, reapplying filter/alignment", flush=True)
        blast_df = pd.read_csv(raw_tsv, sep="\t")
        blast_df["bitscore"] = pd.to_numeric(blast_df["bitscore"], errors="coerce")

    blast_filtered = blast_df[blast_df["bitscore"] >= bitscore_threshold]
    blast_filtered.to_csv(protein_dir / f"{stem}_blast_filtered.tsv", sep="\t", index=False)

    records: list[tuple[str, str]] = [(query_id, seq)]
    seen_seqs = {seq}
    for sseqid in blast_filtered["sseqid"]:
        s = seq_index.get(sseqid)
        if s and s not in seen_seqs:
            records.append((sseqid, s))
            seen_seqs.add(s)

    n_hits = len(records) - 1
    print(f"    {len(blast_df)} raw hits → {n_hits} after bitscore≥{bitscore_threshold}", flush=True)

    if n_hits == 0:
        print(f"    no hits after filtering — skipping alignment")
        return

    unaligned_fasta = protein_dir / f"{stem}_unaligned.fasta"
    aligned_fasta   = protein_dir / f"{stem}_prophage.fasta"
    _write_fasta(records, unaligned_fasta)
    print(f"    aligning {len(records)} seqs …", flush=True)
    _align_fasta(unaligned_fasta, aligned_fasta)
    n_aligned = sum(1 for _ in SeqIO.parse(aligned_fasta, "fasta"))
    print(f"    saved → {aligned_fasta.name} ({n_aligned} seqs)")


def blast_repr_vs_prophage(
    best_predictors_csv: Path,
    prophage_faa_pattern: str,
    per_locus_dir: Path,
    experimental_fastas: list[Path],
    blast_db_dir: Path,
    run_blast: bool = True,
) -> None:
    """
    BLAST all proteins (predictors + experimental) against prophage proteins.
    Outputs written directly to each protein's folder:

    Predictors (n=18):
        per_locus_dir/{KL}/protein/{stem}_blast_raw.tsv
        per_locus_dir/{KL}/protein/{stem}_blast_filtered.tsv
        per_locus_dir/{KL}/protein/{stem}_unaligned.fasta
        per_locus_dir/{KL}/protein/{stem}_prophage.fasta

    Experimental (n=4):
        per_locus_dir/experimental/{subfolder}/{stem}_blast_raw.tsv
        per_locus_dir/experimental/{subfolder}/{stem}_blast_filtered.tsv
        per_locus_dir/experimental/{subfolder}/{stem}_unaligned.fasta
        per_locus_dir/experimental/{subfolder}/{stem}_prophage.fasta
    """
    blast_db_dir.mkdir(parents=True, exist_ok=True)
    bitscore_threshold = 500

    print(f"  Building prophage BLAST DB …")
    db_path, seq_index = _build_blast_db(prophage_faa_pattern, blast_db_dir)
    print(f"  {len(seq_index)} prophage proteins indexed")

    # --- Predictors ---
    best = pd.read_csv(best_predictors_csv)
    for _, row in best.iterrows():
        locus   = row["locus"]
        pc      = row["PC"]
        rep_seq = row.get("representative_sequence")
        if pd.isna(rep_seq) or not rep_seq:
            print(f"  [{locus}] no representative_sequence — skipping")
            continue
        stem = f"{locus}_{pc}"
        _blast_one_protein(
            query_id=stem, stem=stem, seq=rep_seq,
            db_path=db_path, seq_index=seq_index, db_dir=blast_db_dir,
            protein_dir=per_locus_dir / locus / "protein",
            bitscore_threshold=bitscore_threshold,
            run_blast=run_blast,
        )

    # --- Experimental proteins ---
    for protein_id, info, seq in _load_experimental_records(experimental_fastas):
        _blast_one_protein(
            query_id=protein_id, stem=info["stem"], seq=seq,
            db_path=db_path, seq_index=seq_index, db_dir=blast_db_dir,
            protein_dir=per_locus_dir / "experimental" / info["subfolder"],
            bitscore_threshold=bitscore_threshold,
            run_blast=run_blast,
        )


def export_per_protein(
    best_predictors_csv: Path,
    per_locus_dir: Path,
    experimental_fastas: list[Path],
) -> None:
    """
    Populate protein folders with CIF structure and sequence FASTA.
    (Blast outputs are written directly by blast_repr_vs_prophage.)

    Predictors → per_locus_dir/{KL}/protein/
    Experimental → per_locus_dir/{folder}/protein/{protein_id}/
    """
    # --- Predictors ---
    best = pd.read_csv(best_predictors_csv)
    for _, row in best.iterrows():
        locus      = row["locus"]
        pc         = row["PC"]
        protein_id = row.get("representative_protein_id", "")
        rep_seq    = row.get("representative_sequence", "")
        stem       = f"{locus}_{pc}"

        protein_dir = per_locus_dir / locus / "protein"
        protein_dir.mkdir(parents=True, exist_ok=True)

        if pd.notna(rep_seq) and rep_seq:
            header = protein_id if pd.notna(protein_id) and protein_id else stem
            with open(protein_dir / f"{stem}_sequence.fasta", "w") as fh:
                fh.write(f">{header}\n{rep_seq}\n")
        else:
            print(f"  [warn] {locus}: no representative_sequence")

        print(f"  {locus}: protein/ populated")

    # --- Experimental proteins ---
    for protein_id, info, seq in _load_experimental_records(experimental_fastas):
        protein_dir = per_locus_dir / "experimental" / info["subfolder"]
        protein_dir.mkdir(parents=True, exist_ok=True)

        with open(protein_dir / f"{info['stem']}_sequence.fasta", "w") as fh:
            fh.write(f">{protein_id}\n{seq}\n")

        print(f"  {protein_id}: experimental/{info['subfolder']}/ populated")
