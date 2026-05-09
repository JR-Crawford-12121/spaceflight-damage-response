# Torchlight notebook ↔ OSDR files (`Torchlight_Hackathon_2026.ipynb`)

The starter notebook **`data/raw/Torchlight_Hackathon_2026.ipynb`** is the authoritative description of how each NASA OSD processed Excel file should be read (`pd.read_excel` kwargs, contrasts, and biology notes).

This repository’s **DNA damage-response RNA pathway scorer** intentionally loads **only two** gene-expression tables that match the notebook’s RNA-seq sections:

| Notebook section (conceptual) | Study | Processed file (`file=` in OSDR URL) | Official contrast (from notebook text) | Loaded by `python main.py`? |
|-------------------------------|-------|--------------------------------------|------------------------------------------|------------------------------|
| Whole blood — total RNA-seq gene expression | OSD-569 | `GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx` | **(R+82) vs (L-92, L-44, L-3)** — DESeq2 / pipeline-transcriptome-de | **Yes** (`src/load_data.py`) |
| Whole blood — m6A modification | OSD-569 | `GLDS-561_directm6Aseq_Direct_RNA_seq_m6A_Processed_Data.xlsx` | **(R+1) vs (L-92, L-44, L-3)** — transcript × methylation probability | **No** (different modality; transcript-centric index `transcript_ENSEMBL`) |
| PBMC — snRNA-seq gene expression | OSD-570 | `GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx` | **(R+45) vs (R+1)** — Seurat FindMarkers per cell type | **Yes** (`src/load_data.py`) |
| PBMC — snATAC-seq | OSD-570 | `GLDS-562_snATAC-Seq_PBMC_Chromatin_Accessibility_snATAC-seq_Processed_Data.xlsx` | **(R+1) vs (L-92, L-44, L-3)** — chromatin peaks (coordinates), not gene symbols | **No** (not a gene-row DE table for HGNC pathway overlap) |
| PBMC — V(D)J | OSD-570 | `GLDS-562_scRNA-Seq_VDJ_Results.xlsx` | Repertoire / clonotypes | **No** |
| Spatial transcriptomics (skin) | OSD-570 | `GLDS-566_SpatialTranscriptomics_Skin_Biopsy_Spatially_Resolved_Transcriptomics_Processed_Data.xlsx` | (see notebook) | **No** |

Additional OSD studies (e.g. OSD-571 metabolomics / proteomics) appear later in the notebook; they are **out of scope** for this RNA pathway prototype unless explicitly added later.

## Exact `read_excel` parameters (from notebook code cells)

Parameters are duplicated in **`src/notebook_manifest.py`** as `READ_EXCEL_BY_FILE` for quick parity checks.

### OSD-569 — RNA gene expression (implemented)

```python
rna_blood = pd.read_excel(
    "...GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx",
    skiprows=[0, 1, 2, 3, 4, 5, 6, 9],
    header=[0, 1],
    index_col=0,
)
rna_blood.columns = [f"{str(a)}_{str(b)}".strip("_") for a, b in rna_blood.columns]
```

### OSD-569 — m6A (documented only)

```python
m6A = pd.read_excel(
    "...GLDS-561_directm6Aseq_Direct_RNA_seq_m6A_Processed_Data.xlsx",
    skiprows=[0, 1, 2, 3, 4, 5, 6, 7],
    header=[0, 1],
    index_col=0,
)
m6A.columns = [f"{str(a)}_{str(b)}".strip("_") for a, b in m6A.columns]
m6A.index.name = "transcript_ENSEMBL"
```

### OSD-570 — PBMC snRNA-seq (implemented)

```python
snrnaseq = pd.read_excel(
    "...GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx",
    skiprows=[0, 1, 2, 3, 4, 5, 6],
    index_col=0,
)
```

### OSD-570 — PBMC snATAC-seq (documented only)

```python
snatacseq = pd.read_excel(
    "...GLDS-562_snATAC-Seq_PBMC_Chromatin_Accessibility_snATAC-seq_Processed_Data.xlsx",
    skiprows=[0, 1, 2, 3, 4, 5, 6],
    index_col=0,
)
```

## Why “replanning” does not mean the RNA loaders were wrong

The **gene-level RNA-seq** and **PBMC snRNA-seq** loaders in **`src/load_data.py`** use the **same** `skiprows`, `header`, and `index_col` as the notebook snippets above.  

What was missing was an explicit **inventory** of *other* notebook tables (m6A, ATAC, …) so it is obvious they are **different assays** with **different row semantics**, not drop-in replacements for the pathway RNA workflow.

## Related files

- **`src/notebook_manifest.py`** — machine-readable kwargs + status strings  
- **`references/PROVENANCE.md`** — short pointer + OSDR URL pattern  
- **`docs/ENSEMBL_MAPPING.md`** — Ensembl→HGNC mapping for OSD-569 RNA index only  
