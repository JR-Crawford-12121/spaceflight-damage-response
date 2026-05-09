"""
Inventory of `data/raw/Torchlight_Hackathon_2026.ipynb` OSDR Excel loaders.

Source of truth for `pd.read_excel(...)` kwargs is the notebook JSON (grep/read_excel).
This module documents **parity** with the notebook and which assets `main.py` actually loads.

The DNA damage-response RNA pathway pipeline uses **only**:
  - OSD-569: whole-blood **gene-level RNA-seq** DE table
  - OSD-570: PBMC **snRNA-seq** gene-level DE table

Other tables (m6A, ATAC, VDJ, spatial, OSD-571, …) are listed for planning only.
"""

from __future__ import annotations

from typing import Any

# Filenames must match notebook URL `file=` parameters exactly.

READ_EXCEL_BY_FILE: dict[str, dict[str, Any]] = {
    # --- Implemented in load_data.py (pathway RNA pipeline) ---
    "GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx": {
        "skiprows": [0, 1, 2, 3, 4, 5, 6, 9],
        "header": [0, 1],
        "index_col": 0,
    },
    "GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx": {
        "skiprows": [0, 1, 2, 3, 4, 5, 6],
        "index_col": 0,
    },
    # --- Notebook-only (same studies; not loaded by main.py v1) ---
    "GLDS-561_directm6Aseq_Direct_RNA_seq_m6A_Processed_Data.xlsx": {
        "skiprows": [0, 1, 2, 3, 4, 5, 6, 7],
        "header": [0, 1],
        "index_col": 0,
    },
    "GLDS-562_snATAC-Seq_PBMC_Chromatin_Accessibility_snATAC-seq_Processed_Data.xlsx": {
        "skiprows": [0, 1, 2, 3, 4, 5, 6],
        "index_col": 0,
    },
    "GLDS-562_scRNA-Seq_VDJ_Results.xlsx": {
        "skiprows": [0, 1, 2],
        "index_col": 0,
    },
    "GLDS-566_SpatialTranscriptomics_Skin_Biopsy_Spatially_Resolved_Transcriptomics_Processed_Data.xlsx": {
        "skiprows": [0, 1, 2, 3, 4, 5, 6],
        "index_col": 0,
    },
}

# Files loaded by load_all_tables() today
PIPELINE_GENE_EXPRESSION_FILES: frozenset[str] = frozenset(
    {
        "GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx",
        "GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx",
    }
)

NOTEBOOK_SUMMARY: list[dict[str, str]] = [
    {
        "study": "OSD-569",
        "modality": "Whole blood RNA-seq — gene expression (DESeq2 / pipeline-transcriptome-de)",
        "file": "GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx",
        "comparison": "(R+82) vs (L-92, L-44, L-3)",
        "pipeline": "loaded — pathway RNA scorer",
    },
    {
        "study": "OSD-569",
        "modality": "m6A modification (transcript-level; m6Anet / methyl-kit)",
        "file": "GLDS-561_directm6Aseq_Direct_RNA_seq_m6A_Processed_Data.xlsx",
        "comparison": "(R+1) vs (L-92, L-44, L-3)",
        "pipeline": "not loaded — different biology & row keys (transcripts)",
    },
    {
        "study": "OSD-570",
        "modality": "PBMC snRNA-seq — gene expression (Seurat FindMarkers)",
        "file": "GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx",
        "comparison": "(R+45) vs (R+1)",
        "pipeline": "loaded — pathway RNA scorer",
    },
    {
        "study": "OSD-570",
        "modality": "PBMC snATAC-seq — chromatin peaks (coordinates)",
        "file": "GLDS-562_snATAC-Seq_PBMC_Chromatin_Accessibility_snATAC-seq_Processed_Data.xlsx",
        "comparison": "(R+1) vs (L-92, L-44, L-3)",
        "pipeline": "not loaded — peaks not HGNC gene rows",
    },
    {
        "study": "OSD-570",
        "modality": "V(D)J repertoire",
        "file": "GLDS-562_scRNA-Seq_VDJ_Results.xlsx",
        "comparison": "(see notebook)",
        "pipeline": "not loaded — clonotypes, not bulk DE genes",
    },
    {
        "study": "OSD-570",
        "modality": "Spatial transcriptomics (skin biopsy)",
        "file": "GLDS-566_SpatialTranscriptomics_Skin_Biopsy_Spatially_Resolved_Transcriptomics_Processed_Data.xlsx",
        "comparison": "(see notebook)",
        "pipeline": "not loaded — separate assay / future plugin",
    },
]


def assert_rnaseq_loader_matches_notebook() -> None:
    """Optional sanity check (e.g. unit test): dict entries match notebook cells."""
    assert READ_EXCEL_BY_FILE[
        "GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx"
    ] == {
        "skiprows": [0, 1, 2, 3, 4, 5, 6, 9],
        "header": [0, 1],
        "index_col": 0,
    }
    assert READ_EXCEL_BY_FILE[
        "GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx"
    ] == {"skiprows": [0, 1, 2, 3, 4, 5, 6], "index_col": 0}
