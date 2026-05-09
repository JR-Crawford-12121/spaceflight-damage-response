"""
Deep inspection of the GLDS-561 multi-sheet processed workbook (OSD-569).

Writes: outputs/logs/osd569_workbook_inspection.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from ensembl_mapping_sources import (  # noqa: E402
    extract_mapping_from_wide_table,
    find_symbol_column,
)
from gene_id_mapping import strip_ensembl_version  # noqa: E402
from load_data import OSD569_FILENAME, locate_file  # noqa: E402
from osd569_io import read_osd569_processed_wide  # noqa: E402
from osd569_parse import (  # noqa: E402
    DEFAULT_OSD569_SHEET,
    extract_comparison_from_sheet_preview,
    find_two_row_header_indices,
    list_sheet_comparisons,
    ordered_adj_p_candidates,
    ordered_effect_candidates,
    ordered_raw_p_candidates,
    read_sheet_two_row_header_mi,
    select_first_numeric_column,
    flatten_multiindex_columns,
)


def _first_n_nonnull(series: pd.Series, n: int = 5) -> list[str]:
    s = pd.to_numeric(series, errors="coerce")
    out: list[str] = []
    for v in s:
        if pd.notna(v):
            out.append(str(v))
        if len(out) >= n:
            break
    return out


def _nonnull_count(series: pd.Series) -> int:
    return int(pd.to_numeric(series, errors="coerce").notna().sum())


def inspect_one_sheet(
    xlsx_path: Path,
    sheet_name: str,
    lines: list[str],
) -> None:
    lines.append("=" * 80)
    lines.append(f"Sheet: {sheet_name!r}")
    lines.append("=" * 80)
    try:
        preview = pd.read_excel(
            xlsx_path, sheet_name=sheet_name, header=None, nrows=50, engine="openpyxl"
        )
    except Exception as exc:
        lines.append(f"  ERROR reading sheet: {exc}")
        lines.append("")
        return

    comp = extract_comparison_from_sheet_preview(preview)
    lines.append(f"comparison (from metadata / top rows, best-effort): {comp!r}")
    h0 = find_two_row_header_indices(preview)
    lines.append(f"detected data header start row (0-based, first of two header rows): {h0!r}")
    if h0 is None:
        lines.append("  Could not auto-detect two-row header; skipping column scan for this sheet.")
        lines.append("")
        return

    try:
        wide_mi = read_sheet_two_row_header_mi(xlsx_path, sheet_name, h0)
        wide, _ = flatten_multiindex_columns(wide_mi)
    except Exception as exc:
        lines.append(f"  ERROR building wide table: {exc}")
        lines.append("")
        return

    cols = list(wide.columns)
    eff_c, _ = select_first_numeric_column(wide, ordered_effect_candidates(cols))
    raw_c, _ = select_first_numeric_column(wide, ordered_raw_p_candidates(cols))
    adj_c, _ = select_first_numeric_column(wide, ordered_adj_p_candidates(cols))
    lines.append("selected columns (numeric-validated priority):")
    lines.append(f"  log2FC / effect: {eff_c!r}")
    lines.append(f"  raw p:         {raw_c!r}")
    lines.append(f"  adj p:         {adj_c!r}")
    lines.append("")

    idx = wide.index.astype(str)[:5].tolist()
    lines.append("first 5 index / Ensembl-style IDs:")
    for v in idx:
        lines.append(f"  {v}")
    lines.append("  (stripped) " + ", ".join(strip_ensembl_version(str(v)) for v in idx))
    lines.append("")

    if eff_c and eff_c in wide.columns:
        lines.append("first 5 numeric log2FC / effect_size values (non-null scan):")
        for v in _first_n_nonnull(wide[eff_c], 5):
            lines.append(f"  {v}")
    else:
        lines.append("first 5 log2FC values: (column missing)")

    if adj_c and adj_c in wide.columns:
        lines.append("first 5 adjusted p-values (non-null scan):")
        for v in _first_n_nonnull(wide[adj_c], 5):
            lines.append(f"  {v}")
    else:
        lines.append("first 5 adjusted p-values: (column missing)")

    lines.append("non-null counts (numeric coerce) for priority-selected columns:")
    for name, c in (("effect", eff_c), ("raw_p", raw_c), ("adj_p", adj_c)):
        if c and c in wide.columns:
            lines.append(f"  {name} ({c}): {_nonnull_count(wide[c])} / {len(wide)}")
    lines.append("")

    for c in cols:
        s = wide[c]
        if pd.api.types.is_numeric_dtype(s):
            continue
        coerced = pd.to_numeric(s, errors="coerce")
        nn = int(coerced.notna().sum())
        if nn == 0:
            continue
        nlow = c.lower()
        if "log2" in nlow or "p" in nlow or "fc" in nlow:
            lines.append(f"  candidate numeric column {c!r}: non-null {nn} / {len(wide)}")

    lines.append("")


def main() -> int:
    out_path = ROOT / "outputs" / "logs" / "osd569_workbook_inspection.txt"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [
        "OSD-569 GLDS-561 processed workbook — full inspection",
        f"repo root: {ROOT}",
        "",
    ]

    xlsx_path = locate_file(OSD569_FILENAME, ROOT)
    if xlsx_path is None:
        msg = f"ERROR: {OSD569_FILENAME} not found under data/raw/ or data/processed/."
        out_path.write_text("\n".join([msg]), encoding="utf-8")
        print(msg)
        return 1

    lines.append(f"Workbook: {xlsx_path}")
    lines.append("")

    try:
        xl = pd.ExcelFile(xlsx_path, engine="openpyxl")
        sheet_names = xl.sheet_names
    except Exception as exc:
        lines.append(f"ERROR: {exc}")
        out_path.write_text("\n".join(lines), encoding="utf-8")
        return 1

    lines.append("All sheet names:")
    for sn in sheet_names:
        lines.append(f"  - {sn!r}")
    lines.append("")

    cm = list_sheet_comparisons(xlsx_path, sheet_names)
    lines.append("Best-effort comparison string per sheet (metadata scan):")
    for sn in sheet_names:
        lines.append(f"  {sn!r}: {cm.get(sn)!r}")
    lines.append("")

    for sn in sheet_names:
        inspect_one_sheet(xlsx_path, sn, lines)

    lines.append("=" * 80)
    lines.append(f"Default pipeline sheet ({DEFAULT_OSD569_SHEET!r}) — same loader as main.py")
    lines.append("=" * 80)
    wb_map: dict[str, str] = {}
    try:
        wide_default = read_osd569_processed_wide(xlsx_path)
        lines.append(f"shape {wide_default.shape}; columns (first 30): {list(wide_default.columns)[:30]}")
        sym_col = find_symbol_column(wide_default)
        lines.append(f"heuristic symbol column: {sym_col!r}")
        wb_map = extract_mapping_from_wide_table(wide_default)
        lines.append(f"mapping rows extracted from workbook text columns: {len(wb_map)}")
    except Exception as exc:
        lines.append(f"ERROR reading default sheet wide table: {exc}")

    if wb_map and len(wb_map) >= 10:
        proc_dir = ROOT / "data" / "processed"
        proc_dir.mkdir(parents=True, exist_ok=True)
        map_path = proc_dir / "ensembl_to_symbol.tsv"
        rows = [{"ensembl_gene_id": k, "hgnc_symbol": v} for k, v in sorted(wb_map.items())]
        pd.DataFrame(rows).to_csv(map_path, sep="\t", index=False)
        lines.append("")
        lines.append(f"Wrote optional workbook-derived mapping -> {map_path} ({len(rows)} pairs)")
    elif wb_map:
        lines.append("")
        lines.append(
            f"NOTE: Only {len(wb_map)} workbook-derived mappings (<10); "
            "skipped auto-writing ensembl_to_symbol.tsv."
        )

    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote: outputs/logs/osd569_workbook_inspection.txt")
    print(f"Full path: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
