"""
OSD-569 GLDS-561 Excel: per-sheet layout detection, header rows, and DE column selection.

Uses numeric validation: a column is accepted only if pd.to_numeric yields non-null values.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

# Intended contrast for whole-blood RNA Gene Expression processed table (notebook text).
EXPECTED_COMPARISON_CANONICAL = "(R+82)_vs_(L-92_L-44_L-3)"
DEFAULT_OSD569_SHEET = "I4-LP2"

OSD569_DATASET_NAME = "OSD-569_whole_blood_RNA"


def normalize_column_name(name: str) -> str:
    """Lowercase, spaces → underscores, collapse repeated underscores."""
    s = str(name).strip().lower()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_")


def normalize_comparison_text(s: str) -> str:
    """Loose normalization for comparing contrast strings from metadata vs canonical label."""
    t = str(s).lower()
    t = t.replace(",", " ")
    t = re.sub(r"\s+", "", t)
    t = t.replace("(+", "_plus_").replace("(-", "_minus_")
    t = t.replace("(", "").replace(")", "")
    t = t.replace("+", "plus").replace("-", "minus")
    return t


def comparison_matches_expected(metadata_fragment: str) -> bool:
    if not metadata_fragment or not str(metadata_fragment).strip():
        return False
    m = normalize_comparison_text(metadata_fragment)
    exp = normalize_comparison_text(EXPECTED_COMPARISON_CANONICAL)
    keys = ["r", "82", "l", "92", "44", "3"]
    if all(k in m for k in keys) and all(k in exp for k in keys):
        return True
    return exp[:20] in m or m[:30] in exp


def _row_text_blob(row: pd.Series) -> str:
    parts: list[str] = []
    for v in row.values:
        if pd.isna(v):
            continue
        parts.append(str(v))
    return " ".join(parts)


def extract_comparison_from_sheet_preview(preview: pd.DataFrame, max_rows: int = 20) -> str | None:
    """Scan top rows for 'Comparison:' or (R+82) vs (...) style text."""
    pat = re.compile(
        r"\(R\+82\)\s+vs\s+\([^)]+\)",
        re.IGNORECASE,
    )
    for i in range(min(max_rows, len(preview))):
        blob = _row_text_blob(preview.iloc[i])
        m = pat.search(blob)
        if m:
            return m.group(0).strip()
        if "comparison" in blob.lower() and "vs" in blob.lower():
            low = blob.lower()
            idx = low.find("comparison")
            tail = blob[idx:]
            if "vs" in tail:
                return tail.strip()
    for i in range(min(max_rows, len(preview))):
        blob = _row_text_blob(preview.iloc[i])
        if "r+82" in blob.lower() and "vs" in blob.lower() and "l-92" in blob.lower():
            return blob.strip()[:200]
    return None


def _combined_header_labels(row_a: pd.Series, row_b: pd.Series) -> list[str]:
    labels: list[str] = []
    for j in range(min(len(row_a), len(row_b))):
        a = "" if pd.isna(row_a.iloc[j]) else str(row_a.iloc[j]).strip()
        b = "" if pd.isna(row_b.iloc[j]) else str(row_b.iloc[j]).strip()
        combined = f"{a}_{b}".strip("_")
        labels.append(combined)
    return labels


def _header_pair_scores(labels_norm: list[str]) -> tuple[bool, int]:
    joined = " ".join(labels_norm)
    score = 0
    if any("log2" in x and ("fc" in x or "fold" in x) for x in labels_norm):
        score += 3
    if any("deseq2" in x for x in labels_norm):
        score += 2
    if any("pipeline" in x and "transcriptome" in x for x in labels_norm):
        score += 2
    if "p-value" in joined or "p_value" in joined or "padj" in joined:
        score += 2
    if any("adjusted" in x and "p" in x for x in labels_norm):
        score += 2
    looks = score >= 4
    return looks, score


def find_two_row_header_indices(preview: pd.DataFrame, max_start: int = 35) -> int | None:
    """Return best row index i for two-row MultiIndex header."""
    if len(preview) < 3:
        return None
    n = min(max_start, len(preview) - 2)
    best_i: int | None = None
    best_score = -1
    for i in range(n + 1):
        if i + 1 >= len(preview):
            break
        labels = _combined_header_labels(preview.iloc[i], preview.iloc[i + 1])
        norms = [normalize_column_name(c) for c in labels]
        ok, sc = _header_pair_scores(norms)
        if ok and sc > best_score:
            best_score = sc
            best_i = i
    return best_i


def _dedupe_column_names(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for n in names:
        base = n if n else "column"
        k = base
        if k not in seen:
            seen[k] = 0
            out.append(k)
        else:
            seen[k] += 1
            out.append(f"{base}__dup{seen[k]}")
    return out


def flatten_multiindex_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """
    Flatten MultiIndex columns to strings; dedupe. Returns (new_df, raw_repr list).
    """
    if not isinstance(df.columns, pd.MultiIndex):
        names = [str(c) for c in df.columns]
        names = _dedupe_column_names(names)
        out = df.copy()
        out.columns = names
        raw_repr = [repr(c) for c in df.columns]
        return out, raw_repr

    flat: list[str] = []
    raw_repr: list[str] = []
    for col in df.columns:
        raw_repr.append(repr(col))
        if isinstance(col, tuple):
            parts = []
            for level in col:
                parts.append("" if pd.isna(level) else str(level).strip())
            combined = "_".join(p for p in parts if p).strip("_")
            flat.append(combined if combined else "unnamed")
        else:
            flat.append(str(col).strip())
    flat = _dedupe_column_names(flat)
    out = df.copy()
    out.columns = flat
    return out, raw_repr


def read_sheet_two_row_header_mi(
    path: Path,
    sheet_name: str,
    header_start: int,
) -> pd.DataFrame:
    """Read sheet with two-row MultiIndex header (columns not flattened)."""
    return pd.read_excel(
        path,
        sheet_name=sheet_name,
        skiprows=list(range(header_start)),
        header=[0, 1],
        index_col=0,
        engine="openpyxl",
    )


def audit_numeric_column(df: pd.DataFrame, col: str) -> dict[str, Any]:
    """Diagnostics for one column."""
    if col not in df.columns:
        return {
            "column": col,
            "normalized": normalize_column_name(col),
            "non_null_numeric_count": 0,
            "first_5_raw": [],
            "first_5_numeric": [],
            "numeric_min": None,
            "numeric_max": None,
            "error": "column_missing",
        }
    raw = df[col]
    raw_head = [raw.iloc[i] if i < len(raw) else None for i in range(min(5, len(raw)))]
    num = pd.to_numeric(raw, errors="coerce")
    nn = int(num.notna().sum())
    num_head = [num.iloc[i] if i < len(num) else None for i in range(min(10, len(num)))]
    vmin = float(num.min()) if nn else None
    vmax = float(num.max()) if nn else None
    return {
        "column": col,
        "normalized": normalize_column_name(col),
        "non_null_numeric_count": nn,
        "first_5_raw": [repr(x) for x in raw_head],
        "first_10_numeric_converted": [float(x) if pd.notna(x) else None for x in num_head[:10]],
        "numeric_min": vmin,
        "numeric_max": vmax,
    }


def ordered_effect_candidates(columns: list[str]) -> list[tuple[str, str]]:
    """Priority order (column_name, rule_tag)."""
    norms = {c: normalize_column_name(c) for c in columns}
    out: list[tuple[str, str]] = []
    seen: set[str] = set()

    def take(pred: Any, tag: str) -> None:
        for c in columns:
            if c in seen:
                continue
            if pred(norms[c]):
                seen.add(c)
                out.append((c, tag))

    take(lambda n: n == "deseq2_log2fc", "1_DESeq2_log2FC_exact")
    take(lambda n: "deseq2" in n and "log2fc" in n.replace("-", ""), "2_deseq2_and_log2fc")
    take(
        lambda n: "pipeline" in n and "transcriptome" in n and "log2fc" in n.replace("-", ""),
        "3_pipeline_transcriptome_log2FC",
    )
    take(lambda n: "pipeline" in n and "log2fc" in n.replace("-", ""), "4_pipeline_and_log2fc")
    take(
        lambda n: "log2" in n and ("fc" in n or "fold" in n),
        "5_log2_fc_or_fold",
    )
    return out


def ordered_raw_p_candidates(columns: list[str]) -> list[tuple[str, str]]:
    norms = {c: normalize_column_name(c) for c in columns}
    out: list[tuple[str, str]] = []
    seen: set[str] = set()

    def take(pred: Any, tag: str) -> None:
        for c in columns:
            if c in seen:
                continue
            if pred(norms[c]):
                seen.add(c)
                out.append((c, tag))

    take(lambda n: n in ("deseq2_p-value", "deseq2_p_value"), "1_DESeq2_p-value_exact")
    take(
        lambda n: "deseq2" in n
        and ("p-value" in n or "p_value" in n)
        and "adjusted" not in n,
        "2_deseq2_p_not_adjusted",
    )
    take(
        lambda n: "pipeline" in n and "transcriptome" in n and "p-value" in n and "adjusted" not in n,
        "3_pipeline_transcriptome_p",
    )
    take(
        lambda n: "pipeline" in n and ("p-value" in n or "p_value" in n) and "adjusted" not in n,
        "4_pipeline_p_not_adjusted",
    )
    return out


def ordered_adj_p_candidates(columns: list[str]) -> list[tuple[str, str]]:
    norms = {c: normalize_column_name(c) for c in columns}
    out: list[tuple[str, str]] = []
    seen: set[str] = set()

    def take(pred: Any, tag: str) -> None:
        for c in columns:
            if c in seen:
                continue
            if pred(norms[c]):
                seen.add(c)
                out.append((c, tag))

    take(
        lambda n: "deseq2" in n and "adjusted" in n and "p" in n,
        "1_DESeq2_adjusted_p",
    )
    take(lambda n: "deseq2" in n and "adjusted" in n, "2_deseq2_adjusted")
    take(
        lambda n: "pipeline" in n and "transcriptome" in n and "adjusted" in n,
        "3_pipeline_transcriptome_adjusted",
    )
    take(
        lambda n: "pipeline" in n and "adjusted" in n and "p" in n,
        "4_pipeline_adjusted_p",
    )
    return out


def select_first_numeric_column(
    df: pd.DataFrame,
    ordered: list[tuple[str, str]],
) -> tuple[str | None, list[dict[str, Any]]]:
    """
    Try candidates in order; accept first with non-null numeric count > 0.
    Returns (selected_name, audit_rows).
    """
    audits: list[dict[str, Any]] = []
    for col, rule in ordered:
        info = audit_numeric_column(df, col)
        info["priority_rule"] = rule
        nn = info["non_null_numeric_count"]
        info["accepted"] = nn > 0
        info["rejected_reason"] = "" if nn > 0 else "zero_numeric_after_coerce"
        audits.append(info)
        if nn > 0:
            return col, audits
    return None, audits


def multiindex_columns_debug_lines(df: pd.DataFrame) -> list[str]:
    lines: list[str] = []
    lines.append("pandas columns representation:")
    lines.append(f"  type: {type(df.columns).__name__}")
    if isinstance(df.columns, pd.MultiIndex):
        lines.append(f"  MultiIndex names: {df.columns.names}")
        lines.append("  levels[0] sample (first 15): " + repr(list(df.columns.levels[0][:15])))
        if df.columns.nlevels > 1:
            lines.append("  levels[1] sample (first 15): " + repr(list(df.columns.levels[1][:15])))
        lines.append("  first 15 tuples:")
        for i, col in enumerate(df.columns[:15]):
            lines.append(f"    {i}: {repr(col)}")
    else:
        lines.append("  first 15 column names: " + repr(list(df.columns[:15])))
    return lines


def write_osd569_raw_sheet_slice(
    path: Path,
    sheet_name: str,
    header_row_0: int,
    df_mi: pd.DataFrame,
    df_flat: pd.DataFrame,
    candidate_cols: list[str],
) -> None:
    """outputs/logs/osd569_raw_sheet_slice.txt"""
    out_path = Path(__file__).resolve().parent.parent / "outputs" / "logs" / "osd569_raw_sheet_slice.txt"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("OSD-569 raw sheet slice debug")
    lines.append(f"workbook: {path}")
    lines.append(f"sheet: {sheet_name!r}")
    lines.append(f"detected header_row_0 (0-based): {header_row_0}")
    lines.append("")

    preview = pd.read_excel(path, sheet_name=sheet_name, header=None, nrows=16, engine="openpyxl")
    lines.append("--- rows 0–15 with header=None ---")
    lines.append(preview.to_string())
    lines.append("")

    lines.append("--- rows around detected header (preview rows max(0,h0-2) .. h0+12) ---")
    raw_more = pd.read_excel(path, sheet_name=sheet_name, header=None, nrows=min(80, header_row_0 + 25), engine="openpyxl")
    start = max(0, header_row_0 - 2)
    end = min(len(raw_more), header_row_0 + 14)
    lines.append(raw_more.iloc[start:end].to_string())
    lines.append("")

    lines.extend(multiindex_columns_debug_lines(df_mi))
    lines.append("")
    lines.append("--- flattened column names (first 40) ---")
    lines.append(repr(list(df_flat.columns[:40])))
    lines.append("")

    lines.append("--- first 10 data rows for candidate numeric columns ---")
    for c in candidate_cols:
        if c in df_flat.columns:
            sub = df_flat[c].head(10)
            lines.append(f"column {c!r}:")
            lines.append(sub.to_string())
        else:
            lines.append(f"column {c!r}: MISSING")
        lines.append("")

    out_path.write_text("\n".join(lines), encoding="utf-8")


def list_sheet_comparisons(path: Path, sheet_names: list[str]) -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for sn in sheet_names:
        try:
            preview = pd.read_excel(path, sheet_name=sn, header=None, nrows=25, engine="openpyxl")
            out[sn] = extract_comparison_from_sheet_preview(preview)
        except Exception:
            out[sn] = None
    return out


@dataclass
class Osd569SheetLoadMeta:
    sheet_name: str
    comparison_from_metadata: str | None
    comparison_matches_expected: bool
    header_row_0_index: int
    effect_size_column: str | None
    raw_p_value_column: str | None
    adjusted_p_value_column: str | None
    multiindex_column_repr: list[str] = field(default_factory=list)
    flattened_column_names: list[str] = field(default_factory=list)
    candidate_audit: dict[str, list[dict[str, Any]]] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


def load_osd569_wide_from_workbook(
    path: Path,
    sheet_name: str = DEFAULT_OSD569_SHEET,
) -> tuple[pd.DataFrame, Osd569SheetLoadMeta]:
    """
    Load OSD-569 processed sheet: MultiIndex read, flatten, numeric-validated column picks.
    """
    xl = pd.ExcelFile(path, engine="openpyxl")
    if sheet_name not in xl.sheet_names:
        raise ValueError(
            f"Sheet {sheet_name!r} not in workbook. Available: {xl.sheet_names!r}"
        )

    preview = pd.read_excel(path, sheet_name=sheet_name, header=None, nrows=45, engine="openpyxl")
    comp = extract_comparison_from_sheet_preview(preview)
    comp_ok = comparison_matches_expected(comp) if comp else False
    warnings: list[str] = []
    if comp and not comp_ok:
        warnings.append(
            f"Metadata comparison {comp!r} does not match expected "
            f"{EXPECTED_COMPARISON_CANONICAL!r} — still loading sheet {sheet_name!r}."
        )
        all_cmp = list_sheet_comparisons(path, xl.sheet_names)
        warnings.append("Available sheets (detected comparison snippets):")
        for sn, c in sorted(all_cmp.items()):
            warnings.append(f"  - {sn}: {c!r}")

    h0 = find_two_row_header_indices(preview)
    if h0 is None:
        raise ValueError(
            f"Could not locate two-row DE header in sheet {sheet_name!r}. "
            "Inspect outputs/logs/osd569_workbook_inspection.txt."
        )

    wide_mi = read_sheet_two_row_header_mi(path, sheet_name, h0)
    wide_flat, mi_repr = flatten_multiindex_columns(wide_mi)
    cols = list(wide_flat.columns)

    eff_order = ordered_effect_candidates(cols)
    raw_order = ordered_raw_p_candidates(cols)
    adj_order = ordered_adj_p_candidates(cols)

    eff_col, eff_audit = select_first_numeric_column(wide_flat, eff_order)
    raw_col, raw_audit = select_first_numeric_column(wide_flat, raw_order)
    adj_col, adj_audit = select_first_numeric_column(wide_flat, adj_order)

    candidate_audit = {
        "effect_size": eff_audit,
        "raw_p_value": raw_audit,
        "adjusted_p_value": adj_audit,
    }

    if eff_col is None:
        warnings.append(
            "No effect_size column with numeric values matched priority rules; "
            f"column sample: {cols[:30]!r}"
        )
    if adj_col is None:
        warnings.append("No adjusted p-value column with numeric values matched priority rules.")
    if raw_col is None:
        warnings.append("No raw p-value column with numeric values matched priority rules.")

    meta = Osd569SheetLoadMeta(
        sheet_name=sheet_name,
        comparison_from_metadata=comp,
        comparison_matches_expected=comp_ok,
        header_row_0_index=h0,
        effect_size_column=eff_col,
        raw_p_value_column=raw_col,
        adjusted_p_value_column=adj_col,
        multiindex_column_repr=mi_repr,
        flattened_column_names=cols,
        candidate_audit=candidate_audit,
        warnings=warnings,
    )

    cand_keys = list({eff_col, raw_col, adj_col, "DESeq2_log2FC", "pipeline-transcriptome-de_log2FC"})
    cand_keys = [c for c in cand_keys if c]
    try:
        write_osd569_raw_sheet_slice(path, sheet_name, h0, wide_mi, wide_flat, cand_keys)
    except Exception:
        pass

    return wide_flat, meta


def _accepted_candidate_numeric_count(audit_list: list[dict[str, Any]] | None) -> int:
    """Non-null numeric count from the first accepted candidate row in an audit list."""
    if not audit_list:
        return 0
    for row in audit_list:
        if row.get("accepted"):
            return int(row.get("non_null_numeric_count") or 0)
    return 0


def compute_osd569_numeric_status(
    out_df: pd.DataFrame,
    wide: pd.DataFrame,
    meta: Osd569SheetLoadMeta,
    mapping: dict[str, str],
    *,
    mapping_loaded: bool,
) -> dict[str, Any]:
    """Fields for JSON / debug report."""
    from gene_sets import PATHWAY_GENE_SETS

    union_sym: set[str] = set()
    for genes in PATHWAY_GENE_SETS.values():
        union_sym.update(str(g).strip().upper() for g in genes)

    n_total = len(out_df)
    eff_nn = int(out_df["effect_size"].notna().sum())
    raw_nn = int(out_df["raw_p_value"].notna().sum())
    padj_nn = int(out_df["adjusted_p_value"].notna().sum())

    cand_eff = _accepted_candidate_numeric_count(meta.candidate_audit.get("effect_size"))
    cand_raw = _accepted_candidate_numeric_count(meta.candidate_audit.get("raw_p_value"))
    cand_adj = _accepted_candidate_numeric_count(meta.candidate_audit.get("adjusted_p_value"))

    genes_upper = out_df["gene"].astype(str).str.strip().str.upper()
    curated_mask = genes_upper.isin(union_sym)
    has_curated_pathway_overlap = bool(curated_mask.any())
    mapped_curated_eff = int((curated_mask & out_df["effect_size"].notna()).sum())
    mapped_curated_padj = int((curated_mask & out_df["adjusted_p_value"].notna()).sum())

    alignment_parts: list[str] = []
    if cand_eff > 0 and eff_nn == 0:
        alignment_parts.append("effect_size")
    if cand_raw > 0 and raw_nn == 0:
        alignment_parts.append("raw_p_value")
    if cand_adj > 0 and padj_nn == 0:
        alignment_parts.append("adjusted_p_value")

    alignment_warning = ""
    if alignment_parts:
        alignment_warning = (
            "OSD-569 numeric column was selected but lost during standardization. "
            "Check index alignment. "
            f"Affected: {', '.join(alignment_parts)}."
        )
        warnings.warn(alignment_warning, UserWarning, stacklevel=2)

    # Scoreable: mapping + at least one curated-pathway gene symbol in table + numeric effects
    # on curated genes + non-null effect sizes overall.
    scoreable = bool(
        mapping_loaded
        and has_curated_pathway_overlap
        and mapped_curated_eff > 0
        and eff_nn > 0
    )
    partial_note = ""
    if eff_nn > 0 and padj_nn == 0:
        partial_note = (
            "Adjusted p-values unavailable after numeric coercion; significance_weight defaults to 1.0 "
            "for gene_signal (see utils.significance_weight). OSD-569 remains scoreable when effect "
            "sizes are present."
        )

    loader_warnings = list(meta.warnings)
    if alignment_warning:
        loader_warnings.append(alignment_warning)

    return {
        "selected_sheet": meta.sheet_name,
        "comparison": meta.comparison_from_metadata or "",
        "comparison_matches_notebook_expected": meta.comparison_matches_expected,
        "expected_comparison_canonical": EXPECTED_COMPARISON_CANONICAL,
        "header_row_0_index": meta.header_row_0_index,
        "effect_size_source_column": meta.effect_size_column or "",
        "raw_p_value_source_column": meta.raw_p_value_column or "",
        "adjusted_p_value_source_column": meta.adjusted_p_value_column or "",
        "candidate_effect_non_null_count": cand_eff,
        "candidate_raw_p_non_null_count": cand_raw,
        "candidate_adj_p_non_null_count": cand_adj,
        "standardized_effect_non_null_count": eff_nn,
        "standardized_raw_p_non_null_count": raw_nn,
        "standardized_adj_p_non_null_count": padj_nn,
        "effect_size_non_null_count": eff_nn,
        "raw_p_value_non_null_count": raw_nn,
        "adjusted_p_value_non_null_count": padj_nn,
        "total_rows": n_total,
        "mapped_curated_genes_with_effect_size": mapped_curated_eff,
        "mapped_curated_genes_with_adjusted_p_value": mapped_curated_padj,
        "has_curated_pathway_overlap": has_curated_pathway_overlap,
        "mapping_loaded": mapping_loaded,
        "osd569_scoreable": scoreable,
        "osd569_partial_score_note": partial_note,
        "alignment_warning": alignment_warning,
        "loader_warnings": loader_warnings,
        "candidate_column_audit": meta.candidate_audit,
        "multiindex_column_repr_sample": meta.multiindex_column_repr[:25],
    }


def write_osd569_loader_debug(
    path: Path,
    out_df: pd.DataFrame,
    meta: Osd569SheetLoadMeta,
    numeric_status: dict[str, Any],
) -> None:
    lines: list[str] = []
    lines.append("OSD-569 loader debug (GLDS-561 Gene Expression Processed)")
    lines.append("")
    lines.append("multiindex column repr (first 25, stored at load time):")
    for s in (meta.multiindex_column_repr[:25] if meta.multiindex_column_repr else []):
        lines.append(f"  {s}")
    lines.append("")
    lines.append(f"selected_sheet: {numeric_status.get('selected_sheet')}")
    lines.append(f"sheet_metadata_comparison: {numeric_status.get('comparison')!r}")
    lines.append(f"comparison_matches_expected: {numeric_status.get('comparison_matches_notebook_expected')}")
    lines.append(f"detected_header_row_0_index (0-based in sheet): {numeric_status.get('header_row_0_index')}")
    lines.append(f"effect_size_source_column: {numeric_status.get('effect_size_source_column')}")
    lines.append(f"raw_p_value_source_column: {numeric_status.get('raw_p_value_source_column')}")
    lines.append(f"adjusted_p_value_source_column: {numeric_status.get('adjusted_p_value_source_column')}")
    lines.append("")
    lines.append("=" * 72)
    lines.append("1) Candidate numeric column audit (wide table, before standardization)")
    lines.append("=" * 72)
    lines.append(
        f"candidate effect_size non-null count (accepted column): "
        f"{numeric_status.get('candidate_effect_non_null_count')}"
    )
    lines.append(
        f"candidate raw_p_value non-null count (accepted column): "
        f"{numeric_status.get('candidate_raw_p_non_null_count')}"
    )
    lines.append(
        f"candidate adjusted_p_value non-null count (accepted column): "
        f"{numeric_status.get('candidate_adj_p_non_null_count')}"
    )
    lines.append("")
    lines.append("=== Candidate column audits (effect_size) ===")
    for row in numeric_status.get("candidate_column_audit", {}).get("effect_size", []):
        lines.append(
            f"  {row.get('column')!r} | norm={row.get('normalized')!r} | "
            f"rule={row.get('priority_rule')} | nn={row.get('non_null_numeric_count')} | "
            f"accepted={row.get('accepted')} | {row.get('rejected_reason', '')}"
        )
        lines.append(f"    first_5_raw: {row.get('first_5_raw')}")
        lines.append(f"    first_10_numeric: {row.get('first_10_numeric_converted')}")
    lines.append("")
    lines.append("=== Candidate column audits (raw_p_value) ===")
    for row in numeric_status.get("candidate_column_audit", {}).get("raw_p_value", []):
        lines.append(
            f"  {row.get('column')!r} | nn={row.get('non_null_numeric_count')} | accepted={row.get('accepted')}"
        )
        lines.append(f"    first_5_raw: {row.get('first_5_raw')}")
    lines.append("")
    lines.append("=== Candidate column audits (adjusted_p_value) ===")
    for row in numeric_status.get("candidate_column_audit", {}).get("adjusted_p_value", []):
        lines.append(
            f"  {row.get('column')!r} | nn={row.get('non_null_numeric_count')} | accepted={row.get('accepted')}"
        )
        lines.append(f"    first_5_raw: {row.get('first_5_raw')}")
    lines.append("")
    lines.append("=" * 72)
    lines.append("2) Standardized dataframe validation (after positional assignment)")
    lines.append("=" * 72)
    lines.append(f"total_rows: {numeric_status.get('total_rows')}")
    lines.append(
        f"standardized effect_size non-null count: {numeric_status.get('standardized_effect_non_null_count')}"
    )
    lines.append(
        f"standardized raw_p_value non-null count: {numeric_status.get('standardized_raw_p_non_null_count')}"
    )
    lines.append(
        f"standardized adjusted_p_value non-null count: {numeric_status.get('standardized_adj_p_non_null_count')}"
    )
    lines.append(
        f"mapped_curated_genes_with_effect_size: {numeric_status.get('mapped_curated_genes_with_effect_size')}"
    )
    lines.append(
        f"mapped_curated_genes_with_adjusted_p_value: "
        f"{numeric_status.get('mapped_curated_genes_with_adjusted_p_value')}"
    )
    lines.append(f"has_curated_pathway_overlap: {numeric_status.get('has_curated_pathway_overlap')}")
    lines.append(f"mapping_loaded: {numeric_status.get('mapping_loaded')}")
    lines.append(f"osd569_scoreable: {numeric_status.get('osd569_scoreable')}")
    if numeric_status.get("alignment_warning"):
        lines.append(f"alignment_warning: {numeric_status.get('alignment_warning')}")
    if numeric_status.get("osd569_partial_score_note"):
        lines.append(f"note: {numeric_status.get('osd569_partial_score_note')}")
    lines.append("")
    lines.append("first 10 standardized rows (key columns):")
    cols10 = [
        "original_gene_id",
        "ensembl_id",
        "gene",
        "effect_size",
        "raw_p_value",
        "adjusted_p_value",
    ]
    present = [c for c in cols10 if c in out_df.columns]
    head = out_df[present].head(10)
    lines.append(head.to_string())
    lines.append("")
    if numeric_status.get("loader_warnings"):
        lines.append("warnings:")
        for w in numeric_status["loader_warnings"]:
            lines.append(f"  - {w}")
        lines.append("")
    lines.append("first 20 mapped curated genes (gene column used for pathways):")
    from gene_sets import PATHWAY_GENE_SETS

    union_sym: set[str] = set()
    for genes in PATHWAY_GENE_SETS.values():
        union_sym.update(str(g).strip().upper() for g in genes)
    genes_upper = out_df["gene"].astype(str).str.strip().str.upper()
    sub = out_df[genes_upper.isin(union_sym)].copy()
    sub = sub[
        [
            "original_gene_id",
            "gene",
            "effect_size",
            "adjusted_p_value",
        ]
    ].head(20)
    if sub.empty:
        lines.append("  (none — no overlap between mapped genes and curated pathway symbols)")
    else:
        for _, r in sub.iterrows():
            lines.append(
                f"  {r['original_gene_id']!s} | {r['gene']!s} | "
                f"effect={r['effect_size']!s} | padj={r['adjusted_p_value']!s}"
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
