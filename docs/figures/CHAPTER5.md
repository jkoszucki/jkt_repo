# Chapter 5 — Klebsiella putative capsular polysaccharide acetyltransferases are widespread in prophages and K-loci

Chapter 5 characterises acetyltransferases encoded in K-loci and prophages across the *Klebsiella pneumoniae* species complex, integrating three detection approaches: BLASTP, keyword-based search, and FoldSeek structural comparison.

> **Status:** Planned — no `processing/` or `figures/` scripts yet.

---

## Scope

| Analysis | Description |
|----------|-------------|
| BLASTP detection (K-loci) | Query: 3 experimental acetyltransferases (K1, K2, K57); database: Kaptive3 K-loci reference dataset (n=3264) |
| Keyword detection (K-loci) | Search Kaptive3 K-loci for keywords "transferase" + "acetyl/acyl"; manual verification |
| FoldSeek detection (K-loci) | Query: acetyltransferases detected by BLASTP/keyword; database: Kaptive3 K-loci reference dataset; protT5 embeddings; filter e-value ≤ 0.001 and product name containing "transferase" and "acetyl/acyl" |
| Functional annotation | FoldSeek vs AlphaFold/UniProt50 v6; top 1000 hits per protein; TMscore categories: 1.0–0.75 and 0.75–0.50 |
| Comparison | Prophage-encoded vs K-locus-encoded acetyltransferases |

---

## Input data

| File | Description |
|------|-------------|
| `input_dir/cps/cps.xlsx` | Experimental acetyltransferase metadata (sheet: `modifying_enzymes`) |
| Kaptive3 K-loci reference dataset (n=3264) | K-locus reference GenBanks |
| `output_dir/acetyltransferase/acetyltransferase.tsv` | Acetyltransferase table from chapter 2 (prophage-encoded + GWAS + experimental) |

---

## Functional annotation categories

FoldSeek hits against AlphaFold/UniProt50 v6 are grouped into:
- Acetyltransferase
- Acyltransferase
- Hexapeptide
- Trimeric LpxA
- OPgC
- Other

---

## Outputs

(TBD — to be defined when implementation begins)

---

## Figures produced by `figures/chapter5/`

(TBD — planned)
