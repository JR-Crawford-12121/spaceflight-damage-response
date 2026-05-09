"""
Scan GLDS-561 workbook for columns that look like DE statistics and report numeric content.

Writes: outputs/logs/osd569_numeric_column_search.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from load_data import OSD569_FILENAME, locate_file  # noqa: E402
from osd569_parse import (  # noqa: E402
    audit_numeric_column,
    extract_comparison_from_sheet_preview,
    find_two_row_header_indices,
    flatten_multiindex_columns,
    normalize_column_name,
    read_sheet_two_row_header_mi,
)


def column_matches_keyword(norm: str) -> bool:
    keys = (
        "log2",
        "fc",
        "p-value",
        "p_value",
        "adjusted",
        "deseq2",
        "pipeline",
        "padj",
        "pval",
    )
    return any(k.replace("-", "_") in norm.replace("-", "_") for k in keys)


def scan_sheet(xlsx: Path, sheet: str, lines: list[str]) -> None:
    lines.append("=" * 88)
    lines.append(f"SHEET: {sheet!r}")
    lines.append("=" * 88)
    try:
        preview = pd.read_excel(
            xlsx, sheet_name=sheet, header=None, nrows=45, engine="openpyxl"
        )
    except Exception as exc:
        lines.append(f"ERROR: {exc}")
        lines.append("")
        return

    comp = extract_comparison_from_sheet_preview(preview)
    lines.append(f"metadata comparison (best-effort): {comp!r}")
    h0 = find_two_row_header_indices(preview)
    lines.append(f"detected header_row_0: {h0}")
    if h0 is None:
        lines.append("(skip wide read)")
        lines.append("")
        return

    try:
        wide_mi = read_sheet_two_row_header_mi(xlsx, sheet, h0)
        wide, _ = flatten_multiindex_columns(wide_mi)
    except Exception as exc:
        lines.append(f"ERROR wide read: {exc}")
        lines.append("")
        return

    lines.append(f"n_rows={len(wide)} n_cols={len(wide.columns)}")
    lines.append("")
    for col in wide.columns:
        norm = normalize_column_name(col)
        if not column_matches_keyword(norm):
            continue
        info = audit_numeric_column(wide, str(col))
        lines.append(f"--- column {col!r} ---")
        lines.append(f"  normalized: {info.get('normalized')}")
        lines.append(f"  non_null_numeric: {info.get('non_null_numeric_count')}")
        lines.append(f"  numeric_min: {info.get('numeric_min')}  numeric_max: {info.get('numeric_max')}")
        lines.append(f"  first_5_raw: {info.get('first_5_raw')}")
        lines.append(f"  first_10_numeric: {info.get('first_10_numeric_converted')}")
        lines.append("")


def main() -> int:
    out = ROOT / "outputs" / "logs" / "osd569_numeric_column_search.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [
        "OSD-569 numeric column search (all sheets)",
        "",
    ]
    xlsx = locate_file(OSD569_FILENAME, ROOT)
    if xlsx is None:
        msg = f"File not found: {OSD569_FILENAME}"
        out.write_text(msg, encoding="utf-8")
        print(msg)
        return 1
    lines.append(f"workbook: {xlsx}")
    lines.append("")

    xl = pd.ExcelFile(xlsx, engine="openpyxl")
    for sheet in xl.sheet_names:
        scan_sheet(xlsx, sheet, lines)

    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote: outputs/logs/osd569_numeric_column_search.txt")
    print(f"Full path: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
