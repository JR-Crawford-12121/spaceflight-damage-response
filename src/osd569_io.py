"""
Read OSD-569 processed gene-expression XLSX the same way as load_data.load_osd569_from_xlsx.

Used by mapping helper scripts so extraction matches the main pipeline.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from gene_id_mapping import strip_ensembl_version
from osd569_parse import DEFAULT_OSD569_SHEET, load_osd569_wide_from_workbook


def read_osd569_processed_wide(path: Path, sheet_name: str | None = None) -> pd.DataFrame:
    """
    Wide DE table: index = gene IDs (usually Ensembl), columns from flattened MultiIndex.
    Delegates to osd569_parse (default sheet ``I4-LP2``, header auto-detection).
    """
    sn = sheet_name or DEFAULT_OSD569_SHEET
    wide, _meta = load_osd569_wide_from_workbook(path, sheet_name=sn)
    return wide


def unique_ensembl_ids_from_wide(df: pd.DataFrame) -> list[str]:
    """Unique cleaned Ensembl gene IDs from the wide table index."""
    seen: set[str] = set()
    out: list[str] = []
    for ix in df.index.astype(str):
        key = strip_ensembl_version(ix.strip())
        if key.startswith("ENSG") and key not in seen:
            seen.add(key)
            out.append(key)
    return out
