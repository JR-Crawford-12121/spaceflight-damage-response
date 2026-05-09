"""
Ensembl gene ID cleaning and Ensembl → HGNC symbol table loading.

Does not use the network. Mapping is optional via data/processed/ensembl_to_symbol.tsv.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# Accepted column names (case-insensitive column-name match in TSV/CSV header)
ENSEMBL_COL_CANDIDATES = [
    "ensembl_gene_id",
    "ensembl_id",
    "gene_id",
    "Ensembl",
    "ensembl",
]

SYMBOL_COL_CANDIDATES = [
    "hgnc_symbol",
    "gene_symbol",
    "symbol",
    "external_gene_name",
    "Gene",
    "gene",
]


def is_ensembl_gene_id(gene_id: str) -> bool:
    """True if the string looks like an Ensembl stable gene ID (e.g. ENSG...)."""
    if gene_id is None or (isinstance(gene_id, float) and np.isnan(gene_id)):
        return False
    s = str(gene_id).strip()
    return s.startswith("ENS") and len(s) >= 7


def strip_ensembl_version(gene_id: str | None) -> str:
    """
    Remove Ensembl version suffix (e.g. ENSG00000000003.15 -> ENSG00000000003).

    Non-Ensembl inputs are returned unchanged (trimmed). Missing values -> empty string.
    """
    if gene_id is None or (isinstance(gene_id, float) and np.isnan(gene_id)):
        return ""
    s = str(gene_id).strip()
    if not s:
        return ""
    if not is_ensembl_gene_id(s):
        return s
    if "." in s:
        return s.split(".", 1)[0]
    return s


def _norm_col(name: str) -> str:
    return str(name).strip().lower().replace(" ", "_")


def _pick_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower_map = {_norm_col(c): c for c in df.columns}
    for cand in candidates:
        lc = _norm_col(cand)
        if lc in lower_map:
            return lower_map[lc]
    return None


def load_ensembl_symbol_mapping(path: Path) -> tuple[dict[str, str], dict[str, Any]]:
    """
    Load Ensembl -> HGNC mapping from TSV/CSV.

    Returns (mapping dict keyed by stripped Ensembl ID, meta dict with columns used).
    """
    meta: dict[str, Any] = {
        "path": str(path),
        "found": path.is_file(),
        "ensembl_column": None,
        "symbol_column": None,
        "n_rows_file": 0,
        "n_pairs_loaded": 0,
    }
    if not path.is_file():
        return {}, meta

    try:
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    except Exception:
        df = pd.read_csv(path, sep=",", dtype=str, comment="#")

    meta["n_rows_file"] = len(df)
    if df.shape[1] < 2:
        return {}, meta

    ens_col = _pick_column(df, ENSEMBL_COL_CANDIDATES)
    sym_col = _pick_column(df, SYMBOL_COL_CANDIDATES)
    if ens_col is None:
        ens_col = df.columns[0]
    if sym_col is None:
        sym_col = df.columns[1]

    meta["ensembl_column"] = ens_col
    meta["symbol_column"] = sym_col

    out: dict[str, str] = {}
    for _, row in df.iterrows():
        raw_e = row.get(ens_col)
        raw_s = row.get(sym_col)
        if pd.isna(raw_e) or pd.isna(raw_s):
            continue
        key = strip_ensembl_version(str(raw_e).strip())
        sym = str(raw_s).strip()
        if key and sym:
            out[key] = sym

    meta["n_pairs_loaded"] = len(out)
    return out, meta
