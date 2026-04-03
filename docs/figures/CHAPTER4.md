# Chapter 4 — Klebsiella capsular polysaccharide acetylation is a rapidly changing modification

Chapter 4 analyses the diversity of *Klebsiella pneumoniae* CPS repeating unit structures across all 81 known K-types, characterises O-acetylation and O-pyruvylation modifications, and computes pairwise structural similarity using a path-based fingerprint score. The central finding is that O-acetylation is the dominant source of CPS structural heterogeneity — six of eight documented near-identical serotype pairs differ only in O-acetylation state.

---

> **Data preparation** is done by `processing/cps-proc/` — see `docs/processing/CPS-PROC.md` for pipeline detail, inputs, outputs, and library modules.

---

## Figures produced by `figures/chapter4/`

| Figure | Panel | Plot file | Description |
|--------|-------|-----------|-------------|
| Figure 1 | A | `figure1-panelA.png` | Cumulative count of solved CPS structures over time (1972–present) |
| Figure 2 | A | `figure2-panelA.png` | CPS repeating unit structure drawings per K-type |
| Figure 2 | B | `figure2-panelB.png` | Summary table of modifications per K-type |
| Figure 3 | A | `figure3-panelA.png` | Pairwise structural similarity heatmap (J_total) |
| Figure 3 | B | `figure3-panelB.png` | Jaccard similarity score distribution histogram |
| Figure 4 | A | `figure4-panelA.png` | Modification frequency (OAc/OPy) per K-type |
| Figure 4 | B | `figure4-panelB.png` | Molecular weight distribution histogram |
| — | — | `plots/ktypes_structures/` | Per-K-type structure drawings (SVG + PNG) |
| — | — | `plots/cytoscape/` | Similarity network: `edge.tsv` + `node.tsv` for Cytoscape |
