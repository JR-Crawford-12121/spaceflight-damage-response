"""
Extract Ensembl→symbol mappings from OSD-569 wide tables or Ensembl REST API.

REST lookup is only for optional helper scripts — main.py stays offline.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from typing import Any

import pandas as pd

from gene_id_mapping import strip_ensembl_version

# Normalized substring hints for flattened OSD-569 column names
SYMBOL_COLUMN_HINTS = (
    "symbol",
    "gene_symbol",
    "hgnc",
    "gene_name",
    "external_gene",
    "gene_name",
)

STAT_HINTS = ("log2fc", "log2_fc", "p-value", "p_value", "padj", "fc", "statistic")


def _norm_col(name: str) -> str:
    return str(name).strip().lower().replace(" ", "_")


def find_symbol_column(df: pd.DataFrame) -> str | None:
    """Best-effort column that holds HGNC-like symbols (not DE statistics)."""
    candidates: list[tuple[int, str]] = []
    for c in df.columns:
        s = _norm_col(str(c))
        if any(h in s for h in STAT_HINTS):
            continue
        if not any(h in s for h in SYMBOL_COLUMN_HINTS):
            continue
        # Prefer hgnc_symbol / gene_symbol style
        score = 0
        if "hgnc" in s:
            score += 2
        if "symbol" in s:
            score += 2
        if "gene_name" in s or "gene.name" in s:
            score += 1
        candidates.append((score, str(c)))
    if not candidates:
        return None
    candidates.sort(key=lambda x: (-x[0], x[1]))
    return candidates[0][1]


def extract_mapping_from_wide_table(df: pd.DataFrame) -> dict[str, str]:
    """
    Build Ensembl stable ID -> symbol using row index (gene IDs) + symbol column if present.
    Returns {} if no suitable symbol column.
    """
    sym_col = find_symbol_column(df)
    if sym_col is None:
        return {}

    mapping: dict[str, str] = {}
    for i in range(len(df)):
        raw_id = df.index[i]
        ens = strip_ensembl_version(str(raw_id).strip())
        if not ens.startswith("ENSG"):
            continue
        sym = df.iloc[i][sym_col]
        if pd.isna(sym):
            continue
        s = str(sym).strip()
        if not s or s.startswith("ENS"):
            continue
        mapping[ens] = s
    return mapping


ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id"

HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
    "User-Agent": (
        "spaceflight-damage-response/1.0 "
        "(https://github.com/; offline-first pipeline mapping helper)"
    ),
}


def lookup_symbols_via_ensembl_rest(
    ensembl_ids: list[str],
    *,
    chunk_size: int = 400,
) -> tuple[dict[str, str], list[str]]:
    """
    Batch lookup via Ensembl REST POST /lookup/id.
    Returns (ensembl_id -> display_name/symbol, list of IDs with no usable symbol).
    """
    mapping: dict[str, str] = {}

    for i in range(0, len(ensembl_ids), chunk_size):
        chunk = ensembl_ids[i : i + chunk_size]
        body = json.dumps({"ids": chunk}).encode("utf-8")
        req = urllib.request.Request(
            ENSEMBL_LOOKUP_URL,
            data=body,
            headers=HEADERS,
            method="POST",
        )
        try:
            with urllib.request.urlopen(req, timeout=180) as resp:
                payload: Any = json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")[:500]
            raise RuntimeError(
                f"Ensembl REST HTTP {exc.code}: {detail}"
            ) from exc
        except urllib.error.URLError as exc:
            raise RuntimeError(f"Ensembl REST network error: {exc}") from exc

        # Response: list of gene objects aligned with request order
        if isinstance(payload, list):
            for j, item in enumerate(payload):
                if not isinstance(item, dict):
                    continue
                if item.get("error"):
                    continue
                rid = item.get("id")
                if rid is None and j < len(chunk):
                    rid = chunk[j]
                if rid is None:
                    continue
                key = strip_ensembl_version(str(rid))
                sym = item.get("display_name") or item.get("display_label")
                if sym:
                    mapping[key] = str(sym).strip()
        elif isinstance(payload, dict):
            # Rare: keyed by id
            for k, item in payload.items():
                if isinstance(item, dict):
                    key = strip_ensembl_version(str(item.get("id", k)))
                    sym = item.get("display_name") or item.get("display_label")
                    if sym:
                        mapping[key] = str(sym).strip()

    unmapped = [e for e in ensembl_ids if e not in mapping]
    return mapping, unmapped
