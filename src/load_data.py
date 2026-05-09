"""
Load differential expression tables from data/raw/ (or data/processed/).

Only **two** Torchlight notebook tables are loaded here — the gene-level RNA-seq exports:

- OSD-569 whole blood: `GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx`
- OSD-570 PBMC snRNA-seq: `GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx`

The notebook also defines m6A, snATAC, VDJ, spatial, etc.; those are **not** drop-in for this
pathway scorer (see `docs/NOTEBOOK_DATA_SOURCES.md`, `src/notebook_manifest.py`).

Loader kwargs match `Torchlight_Hackathon_2026.ipynb` for the two files above.
"""

from __future__ import annotations

import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from column_detection import ColumnDetectionReport, detect_columns
from gene_id_mapping import load_ensembl_symbol_mapping, strip_ensembl_version
from gene_sets import PATHWAY_GENE_SETS
from osd569_parse import (
    DEFAULT_OSD569_SHEET,
    Osd569SheetLoadMeta,
    compute_osd569_numeric_status,
    load_osd569_wide_from_workbook,
    write_osd569_loader_debug,
)
from utils import project_root

# Official OSD filenames from the Torchlight Hackathon starter notebook
OSD569_FILENAME = "GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx"
OSD570_FILENAME = "GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx"

OSD569_COMPARISON = "(R+82)_vs_(L-92_L-44_L-3)"

# Default processed sheet (GeneLab multi-sheet workbook); override via OSD569_SHEET env.
DEFAULT_OSD569_SHEET_NAME = DEFAULT_OSD569_SHEET
OSD570_COMPARISON = "(R+45)_vs_(R+1)"

MAPPING_BASENAME = "ensembl_to_symbol.tsv"


def _search_directories(root: Path) -> list[Path]:
    return [root / "data" / "raw", root / "data" / "processed"]


def locate_file(filename: str, root: Path | None = None) -> Path | None:
    """Return first path where filename exists."""
    root = root or project_root()
    for d in _search_directories(root):
        p = d / filename
        if p.is_file():
            return p
    return None


def load_osd569_from_xlsx(
    path: Path,
    mapping: dict[str, str],
    sheet_name: str | None = None,
    *,
    mapping_loaded: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Load official GLDS-561 Gene Expression Processed workbook.

    Uses per-sheet metadata + two-row header detection (see osd569_parse.py).
    Default sheet name: ``I4-LP2`` (Torchlight / GeneLab multi-sheet layout).
    Override with ``sheet_name=`` or environment variable ``OSD569_SHEET``.
    """
    import os

    sn = sheet_name or os.environ.get("OSD569_SHEET") or DEFAULT_OSD569_SHEET_NAME
    raw, meta = load_osd569_wide_from_workbook(path, sheet_name=sn)
    for w in meta.warnings:
        if not w.startswith("_rule"):
            print(f"OSD-569 loader: {w}")
    out, numeric_status = osd569_wide_to_standard(
        raw,
        mapping,
        source_path=str(path),
        sheet_meta=meta,
        mapping_loaded=mapping_loaded,
    )
    dbg_path = project_root() / "outputs" / "logs" / "osd569_loader_debug.txt"
    write_osd569_loader_debug(dbg_path, out, meta, numeric_status)
    return out, numeric_status


def osd569_wide_to_standard(
    df: pd.DataFrame,
    mapping: dict[str, str],
    source_path: str = "",
    *,
    sheet_meta: Osd569SheetLoadMeta | None = None,
    mapping_loaded: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    One row per gene (Ensembl IDs in index). Effect / p-value columns come from
    priority selection in osd569_parse (DESeq2 vs pipeline-transcriptome-de).

    Adds original_gene_id, stripped Ensembl key, gene_symbol_mapped, and sets gene
    to HGNC symbol when mapping exists else cleaned Ensembl ID (for overlap matching).
    """
    meta = sheet_meta or Osd569SheetLoadMeta(
        sheet_name="unknown",
        comparison_from_metadata=None,
        comparison_matches_expected=False,
        header_row_0_index=-1,
        effect_size_column=None,
        raw_p_value_column=None,
        adjusted_p_value_column=None,
    )

    effect_col = meta.effect_size_column
    padj_col = meta.adjusted_p_value_column
    raw_p_col = meta.raw_p_value_column

    primary_note = "DESeq2"
    eff_for_note = effect_col or ""
    if eff_for_note and "pipeline" in eff_for_note.lower():
        primary_note = "pipeline-transcriptome-de"
    elif eff_for_note and "deseq2" in eff_for_note.lower():
        primary_note = "DESeq2"

    # RangeIndex + positional numeric assignment — assigning Series indexed by Ensembl IDs into a
    # default RangeIndex DataFrame aligns on labels and yields all NaN (silent bug).
    n = len(df)
    ri = pd.RangeIndex(n)
    out = pd.DataFrame(index=ri)

    og = df.index.astype(str).to_numpy()
    out["original_gene_id"] = og
    ens_arr = np.array([strip_ensembl_version(str(x)) for x in og], dtype=object)
    out["ensembl_id"] = ens_arr
    gsm = pd.Series(ens_arr).map(lambda k: mapping.get(k))
    out["gene_symbol_mapped"] = gsm.to_numpy()
    out["gene"] = gsm.fillna(pd.Series(ens_arr)).astype(str).to_numpy()

    if effect_col and effect_col in df.columns:
        out["effect_size"] = pd.to_numeric(df[effect_col], errors="coerce").to_numpy()
    else:
        out["effect_size"] = np.full(n, np.nan, dtype=float)

    if padj_col and padj_col in df.columns:
        out["adjusted_p_value"] = pd.to_numeric(df[padj_col], errors="coerce").to_numpy()
    else:
        out["adjusted_p_value"] = np.full(n, np.nan, dtype=float)

    if raw_p_col and raw_p_col in df.columns:
        out["raw_p_value"] = pd.to_numeric(df[raw_p_col], errors="coerce").to_numpy()
    else:
        out["raw_p_value"] = np.full(n, np.nan, dtype=float)

    cmp_label = (
        meta.comparison_from_metadata.strip()
        if meta.comparison_from_metadata
        else OSD569_COMPARISON
    )
    out["comparison"] = cmp_label if cmp_label else OSD569_COMPARISON
    out["dataset_name"] = "OSD-569_whole_blood_RNA"
    out["cell_type"] = "unknown"
    out["timepoint"] = "unknown"
    out["source_file"] = source_path
    out["primary_pipeline"] = primary_note

    processed_path = project_root() / "data" / "processed" / "osd569_rnaseq_standardized.tsv"
    processed_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(processed_path, sep="\t", index=False)

    numeric_status = compute_osd569_numeric_status(
        out, df, meta, mapping, mapping_loaded=mapping_loaded
    )
    return out, numeric_status


def _union_curated_pathway_symbols() -> set[str]:
    """HGNC symbols used in curated pathway gene sets (uppercase)."""
    s: set[str] = set()
    for genes in PATHWAY_GENE_SETS.values():
        for g in genes:
            s.add(str(g).strip().upper())
    return s


def summarize_osd569_mapping_stats(out_df: pd.DataFrame, mapping: dict[str, str]) -> dict[str, Any]:
    """Counts for debug report and JSON (OSD-569 rows only)."""
    stripped = out_df["ensembl_id"].astype(str)
    uniq_stripped = stripped.unique()
    n_total_unique = len(uniq_stripped)
    mapped_keys = {k for k in uniq_stripped if k in mapping and str(mapping[k]).strip()}
    n_mapped_unique = len(mapped_keys)
    n_unmapped_unique = n_total_unique - n_mapped_unique
    pct = (100.0 * n_mapped_unique / n_total_unique) if n_total_unique else 0.0
    examples: list[str] = []
    for k in list(sorted(mapped_keys))[:20]:
        examples.append(f"{k} -> {mapping[k]}")

    union_sym = _union_curated_pathway_symbols()
    # Before mapping: identifiers are cleaned Ensembl IDs (no HGNC symbols) — overlap with curated lists is typically 0.
    genes_before = {str(x).strip().upper() for x in uniq_stripped}
    n_overlap_before = len(genes_before & union_sym)
    # After mapping: standardized `gene` column (symbols where mapped).
    genes_after = {str(x).strip().upper() for x in out_df["gene"].unique()}
    n_overlap_after = len(genes_after & union_sym)

    return {
        "total_unique_ensembl_ids": int(n_total_unique),
        "mapped_unique_ids": int(n_mapped_unique),
        "unmapped_unique_ids": int(n_unmapped_unique),
        "percent_mapped": round(pct, 2),
        "example_mappings": examples[:20],
        "unique_genes_matching_any_curated_pathway_before_mapping": int(n_overlap_before),
        "unique_genes_matching_any_curated_pathway_after_mapping": int(n_overlap_after),
    }


def load_osd570_from_xlsx(path: Path) -> tuple[pd.DataFrame, ColumnDetectionReport]:
    """Mirror Colab: skiprows=[0,1,2,3,4,5,6], index_col=0."""
    raw = pd.read_excel(
        path,
        skiprows=[0, 1, 2, 3, 4, 5, 6],
        index_col=0,
        engine="openpyxl",
    )
    sn = raw.reset_index()
    raw_detection = detect_columns(sn)

    gene_col = "Gene" if "Gene" in sn.columns else sn.columns[0]
    out = pd.DataFrame()
    raw_g = sn[gene_col].astype(str).str.strip()
    out["original_gene_id"] = raw_g
    out["gene"] = raw_g
    out["gene_symbol_mapped"] = raw_g
    out["effect_size"] = pd.to_numeric(sn["avg_log2FC"], errors="coerce")
    if "p_val_adj" in sn.columns:
        out["adjusted_p_value"] = pd.to_numeric(sn["p_val_adj"], errors="coerce")
    else:
        out["adjusted_p_value"] = np.nan
    if "p_val" in sn.columns:
        out["raw_p_value"] = pd.to_numeric(sn["p_val"], errors="coerce")
    else:
        out["raw_p_value"] = np.nan

    ct_col = "Cell Type" if "Cell Type" in sn.columns else None
    out["cell_type"] = sn[ct_col].astype(str) if ct_col else "unknown"
    out["comparison"] = OSD570_COMPARISON
    out["dataset_name"] = "OSD-570_PBMC_snRNA"
    out["timepoint"] = "unknown"
    out["ensembl_id"] = None
    out["source_file"] = str(path)

    processed_path = project_root() / "data" / "processed" / "osd570_snrna_standardized.tsv"
    out.to_csv(processed_path, sep="\t", index=False)

    return out, raw_detection


def build_mock_differential_table() -> pd.DataFrame:
    """
    Toy DE rows using real HGNC symbols from pathway lists so scores are non-empty.
    """
    rows = [
        {
            "gene": "TP53",
            "original_gene_id": "TP53",
            "gene_symbol_mapped": "TP53",
            "effect_size": 1.2,
            "adjusted_p_value": 0.02,
            "raw_p_value": 0.001,
            "comparison": "mock_post_vs_pre",
            "dataset_name": "MOCK",
            "cell_type": "unknown",
            "timepoint": "unknown",
            "ensembl_id": None,
            "source_file": "internal_mock",
        },
        {
            "gene": "CDKN1A",
            "effect_size": 0.9,
            "adjusted_p_value": 0.03,
            "raw_p_value": 0.01,
            "comparison": "mock_post_vs_pre",
            "dataset_name": "MOCK",
            "cell_type": "unknown",
            "timepoint": "unknown",
            "ensembl_id": None,
            "source_file": "internal_mock",
        },
        {
            "gene": "HMOX1",
            "effect_size": 1.5,
            "adjusted_p_value": 0.001,
            "raw_p_value": 1e-5,
            "comparison": "mock_post_vs_pre",
            "dataset_name": "MOCK",
            "cell_type": "unknown",
            "timepoint": "unknown",
            "ensembl_id": None,
            "source_file": "internal_mock",
        },
        {
            "gene": "BAX",
            "effect_size": -0.8,
            "adjusted_p_value": 0.04,
            "raw_p_value": 0.02,
            "comparison": "mock_post_vs_pre",
            "dataset_name": "MOCK",
            "cell_type": "unknown",
            "timepoint": "unknown",
            "ensembl_id": None,
            "source_file": "internal_mock",
        },
        {
            "gene": "CHEK1",
            "effect_size": 0.6,
            "adjusted_p_value": 0.05,
            "raw_p_value": 0.03,
            "comparison": "mock_post_vs_pre",
            "dataset_name": "MOCK",
            "cell_type": "unknown",
            "timepoint": "unknown",
            "ensembl_id": None,
            "source_file": "internal_mock",
        },
        {
            "gene": "ATM",
            "effect_size": 0.4,
            "adjusted_p_value": 0.15,
            "raw_p_value": 0.08,
            "comparison": "mock_post_vs_pre",
            "dataset_name": "MOCK",
            "cell_type": "unknown",
            "timepoint": "unknown",
            "ensembl_id": None,
            "source_file": "internal_mock",
        },
        # PBMC-style mock with cell types for optional figure
        {
            "gene": "TP53",
            "effect_size": 0.5,
            "adjusted_p_value": 0.02,
            "raw_p_value": 0.01,
            "comparison": "(R+45)_vs_(R+1)",
            "dataset_name": "MOCK_PBMC",
            "cell_type": "B Cell",
            "timepoint": "unknown",
            "ensembl_id": None,
            "source_file": "internal_mock",
        },
        {
            "gene": "CDKN1A",
            "effect_size": 0.7,
            "adjusted_p_value": 0.01,
            "raw_p_value": 0.005,
            "comparison": "(R+45)_vs_(R+1)",
            "dataset_name": "MOCK_PBMC",
            "cell_type": "CD8 T Cell",
            "timepoint": "unknown",
            "ensembl_id": None,
            "source_file": "internal_mock",
        },
    ]
    df = pd.DataFrame(rows)
    for col in ("original_gene_id", "gene_symbol_mapped"):
        if col not in df.columns:
            df[col] = df["gene"]
        else:
            df[col] = df[col].fillna(df["gene"])
    return df


def missing_file_message() -> str:
    return (
        "\nExpected OSD differential expression exports are not found in data/raw/ "
        "or data/processed/.\n\n"
        "Quick option: copy Torchlight_Hackathon_2026.ipynb into data/raw/, then run:\n"
        "  python main.py --download-from-notebook\n"
        "That parses the official OSDR URLs from the notebook and downloads:\n"
        f"  - {OSD569_FILENAME}\n"
        f"  - {OSD570_FILENAME}\n\n"
        "Manual option: download using the Colab / references/PROVENANCE.md and place "
        "those filenames in data/raw/.\n\n"
        "Optional: add data/processed/ensembl_to_symbol.tsv (Ensembl ID column + "
        "HGNC symbol column; see scripts/create_ensembl_mapping_template.py).\n\n"
        "Running mock demonstration data so the pipeline still completes.\n"
    )


def load_all_tables(root: Path | None = None) -> tuple[pd.DataFrame, dict[str, Any], bool]:
    """
    Load OSD-569 / OSD-570 when present; otherwise return mock data.

    Returns (combined_df, load_report, used_mock).
    """
    root = root or project_root()
    report: dict[str, Any] = {
        "osd569_path": None,
        "osd570_path": None,
        "mapping_path": None,
        "mapping_loaded": False,
        "errors": [],
        "notes": [],
    }

    mapping_path = root / "data" / "processed" / MAPPING_BASENAME
    report["mapping_path"] = str(mapping_path)
    report["mapping_file_found"] = mapping_path.is_file()
    mapping: dict[str, str] = {}
    mapping_meta: dict[str, Any] = {}
    if mapping_path.is_file():
        mapping, mapping_meta = load_ensembl_symbol_mapping(mapping_path)
        report["ensembl_mapping_file_meta"] = mapping_meta
        report["mapping_loaded"] = mapping_meta.get("n_pairs_loaded", 0) > 0
        report["notes"].append(
            f"Loaded {mapping_meta.get('n_pairs_loaded', 0)} Ensembl→symbol pairs "
            f"from {MAPPING_BASENAME} (columns {mapping_meta.get('ensembl_column')!s} → "
            f"{mapping_meta.get('symbol_column')!s})."
        )
    else:
        report["mapping_loaded"] = False
        report["ensembl_mapping_file_meta"] = {"found": False}

    mapping_loaded_flag = bool(report.get("mapping_loaded"))
    frames: list[pd.DataFrame] = []

    p569 = locate_file(OSD569_FILENAME, root)
    p570 = locate_file(OSD570_FILENAME, root)

    if p569 is not None and not mapping_path.is_file():
        report["notes"].append(
            "OSD-569 whole-blood RNA uses Ensembl IDs; without "
            f"{MAPPING_BASENAME}, pathway overlap with HGNC gene sets is usually zero."
        )

    if p569 is not None:
        report["osd569_path"] = str(p569)
        try:
            df569, osd569_numeric_status = load_osd569_from_xlsx(
                p569, mapping, mapping_loaded=mapping_loaded_flag
            )
            frames.append(df569)
            report["osd569_numeric_status"] = osd569_numeric_status
            report["osd569_ensembl_mapping"] = summarize_osd569_mapping_stats(
                df569, mapping
            )
        except Exception as exc:
            report["errors"].append(f"OSD-569 load failed: {exc}")
            report["errors"].append(traceback.format_exc())

    if p570 is not None:
        report["osd570_path"] = str(p570)
        try:
            df570, det570 = load_osd570_from_xlsx(p570)
            frames.append(df570)
            report["osd570_raw_column_detection"] = det570.to_dict()
        except Exception as exc:
            report["errors"].append(f"OSD-570 load failed: {exc}")
            report["errors"].append(traceback.format_exc())

    used_mock = False
    if not frames:
        print(missing_file_message())
        frames.append(build_mock_differential_table())
        used_mock = True
        report["notes"].append("Used internal mock differential expression table.")

    combined = pd.concat(frames, ignore_index=True)

    if "original_gene_id" not in combined.columns:
        combined["original_gene_id"] = combined["gene"]
    else:
        combined["original_gene_id"] = combined["original_gene_id"].fillna(
            combined["gene"]
        )

    if "gene_symbol_mapped" not in combined.columns:
        combined["gene_symbol_mapped"] = combined["gene"]

    # Harmonize gene identifiers for overlap (HGNC uppercase or Ensembl if unmapped)
    combined["gene"] = combined["gene"].astype(str).str.strip().str.upper()

    return combined, report, used_mock
