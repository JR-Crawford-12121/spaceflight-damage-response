# Spaceflight DNA Damage Response Classifier

## BioBrief

This project analyzes NASA OSDR Inspiration4 omics outputs for exploratory RNA-expression signatures in DNA damage-response biology. It combines OSD-569 whole-blood RNA-seq and OSD-570 PBMC single-nuclei RNA-seq, standardizes heterogeneous tables, maps Ensembl IDs to HGNC gene symbols when needed, scores curated DNA repair / oxidative stress / p53 / apoptosis / cell-cycle checkpoint gene sets, and exports ranked outputs for follow-up. It does **not** directly detect DNA breaks, DNA lesions, mutations, damaged genes, or clinical injury.

For file-by-file interpretation, see `outputs/README_OUTPUTS.md`.

## Problem Statement

Spaceflight exposes humans to radiation, oxidative stress, immune disruption, and other stressors that may affect genome-maintenance biology. NASA OSDR data provide rich astronaut omics outputs, but those outputs are difficult to compare directly because they use different tissue types, contrasts, gene identifiers, and table schemas.

This project asks a focused question:

> Can RNA-expression data highlight DNA damage-response genes and pathways that deserve closer follow-up?

The answer is not a diagnosis. The tool does not claim that a gene is damaged. Instead, it points to DNA damage-response “neighborhoods” where gene-expression activity changed enough to warrant closer inspection.

## Approach

The pipeline follows a reproducible workflow:

1. Load approved NASA OSDR Inspiration4 datasets.
2. Standardize heterogeneous gene-expression tables into one schema.
3. Map OSD-569 Ensembl gene IDs to HGNC symbols for selected pathway genes.
4. Score curated pathway gene sets:

   * `DNA_repair`
   * `oxidative_stress`
   * `p53_pathway`
   * `apoptosis`
   * `cell_cycle_checkpoint`
5. Separate pathway-level signals from weak gene-level observations.
6. Validate OSD-569 row identity so Ensembl IDs stay matched to the correct gene symbols.
7. Export CSV, JSON, Markdown, and heatmap outputs.

The main ranking statistic is `pathway_rank_metric`. It is an exploratory ranking metric, not a clinical severity score and not a DNA lesion burden.

Evidence levels are separated by overlap:

* `0` overlapping genes: no pathway support.
* `1–2` overlapping genes: gene-level observation only.
* `3+` overlapping genes: possible pathway-level signal.

## Datasets Used

### OSD-569: whole-blood RNA-seq

* Tissue/source: whole blood.
* Data type: direct RNA-seq / differential gene expression.
* Selected sheet: `I4-LP2`.
* Contrast: `R+82` vs pre-flight timepoints `L-92, L-44, L-3`.
* Gene identifier format: Ensembl IDs.
* Required processing: Ensembl-to-HGNC mapping for curated pathway genes.

OSD-569 is now scoreable in this pipeline. Mapping, numeric parsing, and row-identity validation all pass. It contributes broad whole-blood pathway coverage, but its observed pathway effects are small and do not currently pass the main-candidate threshold.

### OSD-570: PBMC snRNA-seq

* Tissue/source: PBMC immune cells.
* Data type: single-nuclei RNA-seq differential expression.
* Contrast used by the processed table: `R+45` vs `R+1`.
* Important caveat: this is recovery-period timing, not direct flight-vs-pre-flight baseline.
* Gene identifier format: HGNC-like gene symbols.

OSD-570 provides the strongest main candidate pathway-level signals in the current run.

## Key Results

The strongest candidate pathway-level signals come from OSD-570 PBMC data, especially:

* `DNA_repair` in the `Other` cell group.
* `oxidative_stress` in `CD16+ Monocyte`.
* `cell_cycle_checkpoint` in `CD16+ Monocyte`.
* `oxidative_stress` in `Dendritic Cell`.
* `DNA_repair` in `CD16+ Monocyte`.

The strongest pathway-level candidate signal is:

* Dataset: `OSD-570_PBMC_snRNA`
* Comparison: `(R+45)_vs_(R+1)`
* Cell type: `Other`
* Pathway: `DNA_repair`
* `pathway_rank_metric`: about `1.619`
* Overlap: `3` curated genes

OSD-569 now contributes scoreable whole-blood pathway rows with stronger gene-set coverage, including broad overlap with curated DNA damage-response gene sets. However, its mean signed effects are small, and no OSD-569 pathway currently passes the main-candidate threshold.

## Validation Checks

The pipeline includes several validation checks to avoid misleading outputs.

### OSD-569 mapping and numeric parsing

* Ensembl-to-HGNC mapping loaded: yes.
* Mapping pairs loaded: `63`.
* Mapped curated OSD-569 genes: `63`.
* OSD-569 scoreable: yes.
* OSD-569 selected sheet: `I4-LP2`.
* OSD-569 effect-size column: `DESeq2_log2FC`.
* OSD-569 adjusted p-value column: `DESeq2_adjusted p-value`.
* OSD-569 effect-size non-null count: `61,852`.
* OSD-569 adjusted p-value non-null count: `18,418`.

### OSD-569 mapping integrity

The pipeline checks that OSD-569 top-gene rows preserve the correct Ensembl ID to HGNC symbol mapping.

Integrity check result:

* Checked OSD-569 top-gene rows: `77`.
* Mismatched rows: `0`.
* Rows removed from export: `0`.
* Passed: yes.

Canonical spot checks passed:

* `ENSG00000012048 → BRCA1`
* `ENSG00000079246 → XRCC5`
* `ENSG00000012061 → ERCC1`
* `ENSG00000141510 → TP53`
* `ENSG00000149311 → ATM`

### Reproducibility

The core CSV, JSON, and narrative outputs matched across two runs from the same inputs.

Matched files:

* `outputs/tables/pathway_scores_evidence_filtered.csv`
* `outputs/tables/main_candidate_signals.csv`
* `outputs/tables/top_genes_from_pathway_level_signals.csv`
* `outputs/tables/top_damage_response_genes.csv`
* `outputs/logs/initial_analysis_summary.json`
* `outputs/logs/analysis_narrative.md`

## Main Outputs

### Judge-facing summary

`outputs/logs/follow_up_priorities.md`

A readable follow-up priorities report (signals worth closer inspection, caveats, and practical framing). Generate after a run with `python scripts/build_follow_up_priorities.py`. For a file-by-file tour of artifacts, see `outputs/README_OUTPUTS.md`.

### Main figure

`outputs/figures/evidence_filtered_pathway_heatmap.png`

A shared-scale heatmap of scoreable pathway rows. Color shows mean signed log2 fold change. Blank cells indicate insufficient overlap. The heatmap is exploratory only and does not show clinical severity.

### Main tables

`outputs/tables/main_candidate_signals.csv`

The most important table for final interpretation. It contains pathway-level candidate signals that pass overlap and ranking filters.

`outputs/tables/pathway_scores_evidence_filtered.csv`

Scoreable pathway rows with sufficient curated-gene overlap.

`outputs/tables/top_genes_from_pathway_level_signals.csv`

Genes supporting pathway-level rows. OSD-569 rows pass Ensembl-to-HGNC integrity validation.

`outputs/tables/top_damage_response_genes.csv`

A broader gene-level ranked table. This includes lower-overlap observations and should not be treated as pathway-level evidence by itself.

`outputs/tables/pathway_scores_overlap_but_unscored.csv`

Rows with overlap but missing scoreability, if any. This is a diagnostic table.

### Logs and diagnostics

`outputs/logs/initial_analysis_summary.json`

Machine-readable summary of datasets, columns, validation status, main signals, and limitations.

`outputs/logs/data_debug_report.txt`

Dataset diagnostics, mapping status, overlap counts, and numeric parsing status.

`outputs/logs/osd569_mapping_integrity_check.txt`

Validation that OSD-569 Ensembl IDs match the exported HGNC gene symbols.

`outputs/logs/osd569_loader_debug.txt`

Detailed debug report for OSD-569 workbook parsing and standardization.

## How to Run

### 1. Set up environment

PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

### 2. Place input files

Place the official OSDR / GeneLab processed XLSX files in:

```text
data/raw/
```

Expected files:

```text
data/raw/GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx
data/raw/GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx
```

The curated OSD-569 mapping file should be present at:

```text
data/processed/ensembl_to_symbol.tsv
```

### 3. Run the main pipeline

```powershell
python main.py
```

### 4. Build follow-up priorities (judge-facing summary)

```powershell
python scripts\build_follow_up_priorities.py
```

### 5. Optional reproducibility check

Run the pipeline twice, save the outputs as `outputs_run1` and `outputs_run2`, then run:

```powershell
python scripts\reproducibility_check.py outputs_run1 outputs_run2
```

## Recommended Reading Order for Judges

1. `README.md`
2. `outputs/logs/follow_up_priorities.md`
3. `outputs/figures/evidence_filtered_pathway_heatmap.png`
4. `outputs/tables/main_candidate_signals.csv`
5. `outputs/logs/initial_analysis_summary.json`
6. `outputs/logs/osd569_mapping_integrity_check.txt`

## What This Project Can and Cannot Claim

### Can claim

* It identifies exploratory RNA-expression signatures in DNA damage-response pathways.
* It highlights genes and pathways worth follow-up.
* It handles heterogeneous gene identifiers across OSDR datasets.
* It separates pathway-level signals from gene-level observations.
* It validates OSD-569 Ensembl-to-HGNC row identity.
* It produces repeatable core outputs from the same inputs.

### Cannot claim

* It does not detect DNA breaks.
* It does not detect DNA lesions.
* It does not identify damaged genes.
* It does not diagnose astronaut injury.
* It does not measure mutation burden.
* It does not prove that RNA-expression changes were caused only by DNA damage.

## Limitations

* RNA expression does not directly measure DNA damage.
* Transcriptional shifts may reflect inflammation, stress response, immune state, or cell-composition changes.
* The Inspiration4 sample size is small, so generalization is limited.
* OSD-570 uses a `R+45 vs R+1` recovery-period contrast, not direct flight-vs-pre-flight baseline.
* OSD-569 currently uses a curated-only Ensembl-to-HGNC map for this project’s selected pathway genes.
* `pathway_rank_metric` is an exploratory ranking statistic, not clinical severity.

## Future Work

* Add direct DNA-damage assays when available.
* Expand from curated-only Ensembl-to-HGNC mapping to genome-wide annotation if broader pathway libraries are added.
* Add clearer flight-vs-baseline contrasts when available.
* Compare RNA-seq signals with ATAC-seq, CBC, cytokine, proteomic, or direct DNA-damage assay layers.
* Add formal pathway enrichment methods after the hackathon.
* Add a lightweight web dashboard for browsing pathways, genes, and validation reports.

## Practical Use Case

The realistic use case is **triage for astronaut health monitoring**, not diagnosis. The tool helps answer: *Which DNA damage-response pathways, genes, or immune-cell contexts deserve closer follow-up after spaceflight?*

In the fire analogy, this project does not find the fire. It identifies neighborhoods where fire-response activity changed. That can help researchers decide where to look next with better tools.

Possible follow-up uses include:

* Prioritizing DNA repair, oxidative-stress, p53, apoptosis, or checkpoint pathways for deeper review.
* Highlighting immune-cell contexts, such as PBMC cell groups, where stress-response signatures appear stronger.
* Helping decide which direct assays should be run next, such as sequencing, DNA lesion assays, comet assay, γH2AX staining, cytogenetics, or protein-level validation.
* Supporting long-term astronaut monitoring by flagging expression patterns that may justify closer follow-up over time.

This project should **not** be used to predict cancer, locate cancer, diagnose tissue injury, or identify damaged genes. Blood and PBMC RNA-expression data can suggest biological response patterns, but they cannot localize damage to an organ or prove future disease risk by themselves.

## Project Claim

This tool narrows attention to DNA damage-response neighborhoods where RNA-expression activity changed. It does not identify damaged genes. It produces validated, reproducible, follow-up candidates for deeper biological investigation.
