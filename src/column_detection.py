"""
Detect differential-expression table columns from heterogeneous OSD/processed outputs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd


GENE_CANDIDATES = [
    "gene",
    "gene_symbol",
    "symbol",
    "Gene",
    "external_gene_name",
    "feature",
    "name",
]

EFFECT_CANDIDATES = [
    "logFC",
    "log2FoldChange",
    "avg_log2FC",
    "log2fc",
    "estimate",
    "effect_size",
    "DESeq2_log2FC",
    "pipeline-transcriptome-de_log2FC",
]

PADJ_CANDIDATES = [
    "padj",
    "p_val_adj",
    "adj.P.Val",
    "FDR",
    "q_value",
    "adjusted_p_value",
    "DESeq2_adjusted p-value",
    "pipeline-transcriptome-de_adjusted p-value",
]

RAW_P_CANDIDATES = [
    "pvalue",
    "P.Value",
    "p_val",
    "pvalue",
    "DESeq2_p-value",
    "pipeline-transcriptome-de_p-value",
]

META_CANDIDATES = {
    "comparison": ["comparison", "contrast"],
    "timepoint": ["timepoint"],
    "sample": ["sample"],
    "subject": ["subject"],
    "crew_member": ["crew_member"],
    "cell_type": ["cell_type", "Cell Type", "celltype", "cluster"],
}


def _norm(name: str) -> str:
    return str(name).strip()


def _find_column(
    columns: list[str],
    candidates: list[str],
) -> str | None:
    """First exact match (case-sensitive first pass, then case-insensitive)."""
    col_set = {_norm(c): c for c in columns}
    for cand in candidates:
        if cand in col_set:
            return col_set[cand]
    lower_map = {_norm(c).lower(): c for c in columns}
    for cand in candidates:
        lc = cand.lower()
        if lc in lower_map:
            return lower_map[lc]
    return None


@dataclass
class ColumnDetectionReport:
    """What we matched and what we did not."""

    gene_column: str | None
    effect_column: str | None
    adjusted_p_column: str | None
    raw_p_column: str | None
    meta_columns: dict[str, str | None] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene_column": self.gene_column,
            "effect_column": self.effect_column,
            "adjusted_p_column": self.adjusted_p_column,
            "raw_p_column": self.raw_p_column,
            "meta_columns": dict(self.meta_columns),
            "notes": list(self.notes),
        }


def detect_columns(df: pd.DataFrame) -> ColumnDetectionReport:
    """
    Inspect dataframe columns (excluding index unless reset_index was called).
    """
    columns = [str(c) for c in df.columns]
    report = ColumnDetectionReport(
        gene_column=_find_column(columns, GENE_CANDIDATES),
        effect_column=_find_column(columns, EFFECT_CANDIDATES),
        adjusted_p_column=_find_column(columns, PADJ_CANDIDATES),
        raw_p_column=_find_column(columns, RAW_P_CANDIDATES),
        meta_columns={},
    )
    for logical, cands in META_CANDIDATES.items():
        report.meta_columns[logical] = _find_column(columns, cands)

    if report.gene_column is None and df.index.name:
        iname = _norm(df.index.name)
        if iname.lower() in {"gene", "symbol"}:
            report.notes.append(f"Gene identifiers appear to be in the index ({iname}).")

    missing = []
    if report.effect_column is None:
        missing.append("effect_size")
    if report.gene_column is None and df.index.name not in ("Gene", "gene"):
        missing.append("gene")
    if missing:
        report.notes.append(
            "Could not auto-detect columns: " + ", ".join(missing) + "."
        )

    return report
