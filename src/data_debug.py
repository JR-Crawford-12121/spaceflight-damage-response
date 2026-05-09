"""
Per-dataset diagnostics: identifier overlap with curated pathways and debug reports.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from gene_sets import PATHWAY_GENE_SETS
from utils import project_root


def _norm(g: str) -> str:
    return str(g).strip().upper()


def union_curated_pathway_symbols() -> set[str]:
    """All HGNC symbols appearing in any starter pathway list."""
    s: set[str] = set()
    for genes in PATHWAY_GENE_SETS.values():
        s.update(_norm(x) for x in genes)
    return s


def dataset_unique_genes(df: pd.DataFrame) -> set[str]:
    return set(df["gene"].map(_norm))


def pathway_set_overlap_counts(dataset_genes: set[str]) -> dict[str, int]:
    """How many genes from the dataset appear in each named pathway list."""
    out: dict[str, int] = {}
    for pname, genes in PATHWAY_GENE_SETS.items():
        gset = {_norm(g) for g in genes}
        out[pname] = len(dataset_genes & gset)
    return out


def count_genes_matching_any_pathway(dataset_genes: set[str]) -> int:
    union_sym = union_curated_pathway_symbols()
    return len(dataset_genes & union_sym)


def datasets_with_zero_pathway_overlap(pathway_scores: pd.DataFrame) -> set[str]:
    """
    Dataset names where every pathway row has zero overlapping genes.

    Such datasets are omitted from pathway plots even if identifiers could match.
    """
    if pathway_scores.empty or "dataset_name" not in pathway_scores.columns:
        return set()
    excluded: set[str] = set()
    for ds in pathway_scores["dataset_name"].astype(str).unique():
        sub = pathway_scores[pathway_scores["dataset_name"].astype(str) == ds]
        mx = pd.to_numeric(sub["overlapping_gene_count"], errors="coerce").max()
        if pd.isna(mx) or float(mx) <= 0:
            excluded.add(str(ds))
    return excluded


def datasets_without_pathway_identifier_overlap(
    standardized_df: pd.DataFrame,
) -> set[str]:
    """
    Dataset names where no gene identifier matches any curated pathway symbol.

    Excludes such datasets from pathway plots (likely Ensembl vs symbol mismatch).
    """
    excluded: set[str] = set()
    union_sym = union_curated_pathway_symbols()
    for ds, grp in standardized_df.groupby("dataset_name"):
        genes = dataset_unique_genes(grp)
        if len(genes & union_sym) == 0:
            excluded.add(str(ds))
    return excluded


def write_data_debug_report(
    standardized_df: pd.DataFrame,
    load_report: dict,
    root: Path | None = None,
    pathway_scores: pd.DataFrame | None = None,
) -> Path:
    """
    Write outputs/logs/data_debug_report.txt with per-input-dataset sections.
    """
    root = root or project_root()
    path = root / "outputs" / "logs" / "data_debug_report.txt"
    path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    lines.append("Spaceflight damage-response pipeline — data debug report")
    lines.append("(standardized columns after loaders; gene column used for overlap = HGNC or cleaned ID)")
    lines.append("")

    mp = load_report.get("mapping_path") if load_report else None
    ml = load_report.get("mapping_loaded") if load_report else False
    mf = load_report.get("mapping_file_found") if load_report else False
    meta = (load_report or {}).get("ensembl_mapping_file_meta") or {}
    n_pairs = int(meta.get("n_pairs_loaded", 0) or 0)

    lines.append("### Global: Ensembl → HGNC mapping file (OSD-569)")
    lines.append(f"  mapping file path: {mp}")
    lines.append(f"  mapping file found on disk: {'yes' if mf else 'no'}")
    lines.append(f"  mapping pairs loaded (dictionary size): {n_pairs}")
    lines.append(f"  mapping usable for OSD-569 (pairs > 0): {'yes' if ml else 'no'}")
    if meta.get("ensembl_column"):
        lines.append(f"  detected Ensembl column: {meta.get('ensembl_column')}")
        lines.append(f"  detected HGNC/symbol column: {meta.get('symbol_column')}")
    if load_report and load_report.get("osd569_ensembl_mapping"):
        om = load_report["osd569_ensembl_mapping"]
        lines.append("  OSD-569 DE table (unique cleaned Ensembl IDs):")
        for key in (
            "total_unique_ensembl_ids",
            "mapped_unique_ids",
            "unmapped_unique_ids",
            "percent_mapped",
        ):
            if key in om:
                lines.append(f"    {key}: {om[key]}")
        if "unique_genes_matching_any_curated_pathway_before_mapping" in om:
            lines.append(
                "    unique genes matching ANY curated pathway symbol "
                "(before mapping / Ensembl-only identifiers): "
                f"{om['unique_genes_matching_any_curated_pathway_before_mapping']}"
            )
        if "unique_genes_matching_any_curated_pathway_after_mapping" in om:
            lines.append(
                "    unique genes matching ANY curated pathway symbol "
                "(after mapping / standardized gene column): "
                f"{om['unique_genes_matching_any_curated_pathway_after_mapping']}"
            )
        if om.get("example_mappings"):
            lines.append("    first 20 example mappings (Ensembl -> symbol):")
            for ex in om["example_mappings"][:20]:
                lines.append(f"      {ex}")

    num569 = (load_report or {}).get("osd569_numeric_status") if load_report else None
    lines.append("")
    lines.append("### OSD-569 numeric status")
    if num569:
        lines.append(f"  selected_sheet: {num569.get('selected_sheet')}")
        lines.append(f"  selected_comparison (metadata): {num569.get('comparison')!r}")
        lines.append(f"  comparison_matches_notebook_expected: {num569.get('comparison_matches_notebook_expected')}")
        lines.append(f"  expected_comparison_canonical: {num569.get('expected_comparison_canonical')}")
        lines.append(f"  effect_size_source_column: {num569.get('effect_size_source_column')}")
        lines.append(f"  adjusted_p_value_source_column: {num569.get('adjusted_p_value_source_column')}")
        lines.append(f"  effect_size_non_null_count: {num569.get('effect_size_non_null_count')}")
        lines.append(f"  adjusted_p_value_non_null_count: {num569.get('adjusted_p_value_non_null_count')}")
        lines.append(
            f"  mapped_curated_genes_with_effect_size: "
            f"{num569.get('mapped_curated_genes_with_effect_size')}"
        )
        lines.append(
            f"  mapped_curated_genes_with_adjusted_p_value: "
            f"{num569.get('mapped_curated_genes_with_adjusted_p_value')}"
        )
        lines.append(f"  osd569_scoreable: {num569.get('osd569_scoreable')}")
        if num569.get("loader_warnings"):
            lines.append("  loader_warnings:")
            for w in num569["loader_warnings"]:
                lines.append(f"    - {w}")
    else:
        lines.append("  (no osd569_numeric_status — OSD-569 file not loaded or loader did not run)")
    lines.append("")
    lines.append(
        "  Full mapping how-to: docs/ENSEMBL_MAPPING.md and "
        "scripts/create_ensembl_mapping_template.py."
    )
    lines.append("")

    if load_report and not ml and load_report.get("osd569_path"):
        lines.append(
            "WARNING: OSD-569 was loaded but no Ensembl→HGNC mapping file was loaded "
            "(or the file had zero usable pairs). Gene identifiers may remain Ensembl IDs; "
            "pathway overlap with HGNC-based curated lists is expected to be zero until "
            "data/processed/ensembl_to_symbol.tsv is added."
        )
        lines.append("")

    if pathway_scores is not None and not pathway_scores.empty:
        lines.append("### Pathway scoring summary (max overlapping genes per dataset)")
        for ds in sorted(pathway_scores["dataset_name"].astype(str).unique()):
            sub = pathway_scores[pathway_scores["dataset_name"].astype(str) == ds]
            mx = pd.to_numeric(sub["overlapping_gene_count"], errors="coerce").max()
            mx_fmt = int(mx) if pd.notna(mx) else "nan"
            lines.append(f"  {ds}: max overlapping_gene_count across pathway rows = {mx_fmt}")
        if load_report and load_report.get("osd569_path"):
            sub569 = pathway_scores[
                pathway_scores["dataset_name"].astype(str).str.contains("OSD-569", na=False)
            ]
            if not sub569.empty:
                mx569 = pd.to_numeric(
                    sub569["overlapping_gene_count"], errors="coerce"
                ).max()
                if pd.isna(mx569) or float(mx569) <= 0:
                    lines.append(
                        "  OSD-569 figures: this dataset is omitted from pathway heatmaps/bar charts "
                        "when max overlapping_gene_count is 0 (whole blood uses Ensembl IDs until "
                        "ensembl_to_symbol.tsv maps them to HGNC symbols)."
                    )
        lines.append("")

    for ds_name in sorted(standardized_df["dataset_name"].astype(str).unique()):
        sub = standardized_df[standardized_df["dataset_name"].astype(str) == ds_name]
        genes = dataset_unique_genes(sub)
        lines.append("=" * 72)
        lines.append(f"dataset_name: {ds_name}")
        lines.append(f"detected gene column (standardized): gene")
        lines.append(f"detected effect-size column: effect_size")
        lines.append(f"detected adjusted p-value column: adjusted_p_value")
        lines.append(f"number of rows: {len(sub)}")
        lines.append(f"number of unique genes: {len(genes)}")
        uniq_sorted = sorted(genes)
        preview = uniq_sorted[:20]
        lines.append(f"first 20 gene names (alphabetical): {', '.join(preview)}")
        lines.append("overlap count with each curated gene set (unique genes):")
        for pname, cnt in sorted(pathway_set_overlap_counts(genes).items()):
            lines.append(f"  {pname}: {cnt}")
        lines.append(
            f"genes matching ANY curated pathway symbol: {count_genes_matching_any_pathway(genes)}"
        )
        lines.append("")

    if load_report:
        lines.append("load_report (brief)")
        lines.append(str(load_report))

    path.write_text("\n".join(lines), encoding="utf-8")
    return path
