"""
Small shared helpers: paths, directory creation, numeric utilities.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from gene_id_mapping import load_ensembl_symbol_mapping, strip_ensembl_version


def project_root() -> Path:
    """Repository root (parent of src/)."""
    return Path(__file__).resolve().parent.parent


def ensure_output_dirs() -> None:
    """Create outputs/tables, outputs/figures, outputs/logs if missing."""
    root = project_root()
    for sub in ("outputs/tables", "outputs/figures", "outputs/logs", "data/processed"):
        (root / sub).mkdir(parents=True, exist_ok=True)


def safe_float(x: Any) -> float | None:
    """Convert to float; return None if not possible."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def significance_weight(adjusted_p_value: float | None) -> float:
    """
    Weight used in scoring: min(-log10(p), 10), or 1.0 if p missing/invalid.
    """
    if adjusted_p_value is None:
        return 1.0
    p = safe_float(adjusted_p_value)
    if p is None or p <= 0 or p > 1:
        return 1.0
    return float(min(-np.log10(p), 10.0))


def read_gene_mapping(path: Path) -> dict[str, str]:
    """Backward-compatible alias: load Ensembl→symbol mapping from TSV/CSV."""
    mapping, _meta = load_ensembl_symbol_mapping(path)
    return mapping


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_json(obj: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, default=str)
