# Spaceflight DNA Damage Response Classifier — Initial Analysis Narrative

**Opening warning:** Exploratory transcriptional signatures only. Does not detect DNA breaks, lesions, or clinical injury.

## What data was analyzed

Datasets in this run: **OSD-569_whole_blood_RNA, OSD-570_PBMC_snRNA**.
Targets are **OSD-569** whole-blood RNA and **OSD-570** PBMC snRNA when official XLSX files are placed in `data/raw/`.

Mock demonstration mode: **no**.

## What worked

- **OSD-570 PBMC snRNA** uses HGNC-like gene symbols in the processed table, so overlap with curated pathway lists is computable without Ensembl mapping.

- **OSD-569:** Gene-symbol mapping and numeric differential-expression parsing succeeded, so OSD-569 contributes pathway-level whole-blood signals where overlap thresholds are met.

OSD-569 now contributes scoreable whole-blood pathway rows with stronger gene-set coverage, but the observed mean signed effects are small and no OSD-569 pathway currently passes the main-candidate threshold. This summary is exploratory and does not describe clinical severity.

## What needs fixing

No current loader blocker. OSD-569 is now scoreable, but its whole-blood effects are small and do not pass the main-candidate threshold.

## Main pathway-level candidate signals

*Rows with **≥3 overlapping curated genes**, **numeric** `pathway_rank_metric` and `average_signed_effect`, AND (**pathway_rank_metric ≥ 0.5** OR **significant_gene_count ≥ 1**). Mapped-but-unscored overlap rows are in `pathway_scores_overlap_but_unscored.csv`. Primary score-based tables: `pathway_scores_evidence_filtered.csv`, `pathway_scores.csv`.*

**OSD-570 PBMC:** the processed contrast is **R+45 vs R+1** — **recovery-period timing**, not flight vs pre-flight baseline.
**pathway_rank_metric** is an exploratory ranking statistic — **not clinical severity** and **not DNA lesion burden**.
Where overlap counts are **3–5 genes**, pathway-level signals are **limited** by coverage.

- **OSD-570_PBMC_snRNA** | (R+45)_vs_(R+1) | Other | **DNA_repair**: exploratory ranking metric 1.619 (limited pathway-level signal); mean signed log2 FC 0.4135.
- **OSD-570_PBMC_snRNA** | (R+45)_vs_(R+1) | CD16+ Monocyte | **oxidative_stress**: exploratory ranking metric 1.136 (limited pathway-level signal); mean signed log2 FC -0.3466.
- **OSD-570_PBMC_snRNA** | (R+45)_vs_(R+1) | CD16+ Monocyte | **cell_cycle_checkpoint**: exploratory ranking metric 1.066 (limited pathway-level signal); mean signed log2 FC 0.4215.
- **OSD-570_PBMC_snRNA** | (R+45)_vs_(R+1) | Dendritic Cell | **oxidative_stress**: exploratory ranking metric 0.649 (limited pathway-level signal); mean signed log2 FC -0.3709.
- **OSD-570_PBMC_snRNA** | (R+45)_vs_(R+1) | CD16+ Monocyte | **DNA_repair**: exploratory ranking metric 0.357 (limited pathway-level signal); mean signed log2 FC 0.1549.

## Gene-level observations

*Rows with **1–2 overlapping genes** are **gene-level observations only** — **not** pathway-level evidence.*

- **OSD-570_PBMC_snRNA** | Other T Cell | apoptosis: n=1 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | CD8+ T Cell | DNA_repair: n=1 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | CD8+ T Cell | cell_cycle_checkpoint: n=1 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | CD4+ T Cell | apoptosis: n=1 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | CD8+ T Cell | apoptosis: n=1 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | Other T Cell | DNA_repair: n=2 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | Other T Cell | cell_cycle_checkpoint: n=2 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | CD4+ T Cell | DNA_repair: n=1 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | CD4+ T Cell | cell_cycle_checkpoint: n=1 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | B Cell | apoptosis: n=1 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | Natural Killer Cell | cell_cycle_checkpoint: n=2 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | Natural Killer Cell | DNA_repair: n=2 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | B Cell | cell_cycle_checkpoint: n=2 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | B Cell | DNA_repair: n=1 — **gene-level observation only**.
- **OSD-570_PBMC_snRNA** | Other | oxidative_stress: n=1 — **gene-level observation only**.

## Top genes from pathway-level signals

Prefer **`outputs/tables/top_genes_from_pathway_level_signals.csv`** for presentation: genes tied to **parent pathways with ≥3 overlapping genes** only (numeric `gene_signal` and `effect_size` required). **`top_damage_response_genes.csv`** includes **all** gene-level rows (including low-overlap observations).

- **XRCC5** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0002
- **XRCC6** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0002
- **GADD45B** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0002
- **ERCC1** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0001
- **PARP1** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0001
- **XPA** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0001
- **MSH2** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0001
- **ATM** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0001
- **TP53** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0001
- **ATR** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0000
- **DDB2** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0000
- **BRCA1** | DNA_repair | OSD-569_whole_blood_RNA | unknown | gene_signal=0.0000

## Limitations

- Exploratory only; **not diagnostic**.
- **RNA expression does not prove DNA lesions.**
- **Small crew sample size** (e.g. Inspiration4 n=4).
- Immune composition and inflammation may confound pathway readouts.
- **OSD-570 R+45 vs R+1** is **recovery timing**, not direct flight perturbation vs baseline.
- **OSD-569** currently uses a curated-only Ensembl→HGNC map for the project's selected pathway genes.
- **OSD-569** passes mapping and numeric-loader validation, but its whole-blood pathway effects are small and do not pass the main-candidate threshold.

## Next steps

1. Expand OSD-569 from curated-only Ensembl→HGNC mapping to a genome-wide annotation map if broader pathway libraries are added.
2. Add or locate **R+1 vs pre-flight** (or similar) contrasts when available.
3. Compare **flight perturbation** vs **recovery shift** explicitly.
4. Expand curated gene sets using established libraries when appropriate.

---

_Machine-readable outputs: `outputs/logs/initial_analysis_summary.json`, `outputs/logs/data_debug_report.txt`._