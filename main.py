"""
Spaceflight DNA Damage Response Classifier — local analysis entry point.

Usage:
    python main.py

This does not diagnose DNA damage; it summarizes exploratory transcriptional
patterns related to DNA damage-response biology from differential expression.
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

import pandas as pd

# Allow imports from src/ without installing a package
_ROOT = Path(__file__).resolve().parent
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from column_detection import ColumnDetectionReport  # noqa: E402
from damage_score import attach_interpretation_labels, score_all_pathways  # noqa: E402
from data_debug import (  # noqa: E402
    datasets_with_zero_pathway_overlap,
    datasets_without_pathway_identifier_overlap,
    write_data_debug_report,
)
from load_data import load_all_tables  # noqa: E402
from narrative import build_analysis_narrative_md  # noqa: E402
from plotting import run_all_plots  # noqa: E402
from summarize_results import (  # noqa: E402
    build_initial_analysis_summary_json,
    build_osd569_score_inclusion_status,
    pathway_level_top_genes_table,
    prepare_pathway_level_top_genes_export,
    save_summary_outputs,
)
from utils import ensure_output_dirs, project_root  # noqa: E402


def _osd569_has_any_pathway_overlap(pathway_scores: pd.DataFrame) -> bool:
    m = pathway_scores["dataset_name"].astype(str).str.contains("OSD-569", na=False)
    if not m.any():
        return False
    sub = pd.to_numeric(
        pathway_scores.loc[m, "overlapping_gene_count"], errors="coerce"
    ).fillna(0)
    return bool((sub > 0).any())


def _osd569_max_pathway_overlap(pathway_scores: pd.DataFrame) -> int:
    m = pathway_scores["dataset_name"].astype(str).str.contains("OSD-569", na=False)
    if not m.any():
        return 0
    v = pd.to_numeric(pathway_scores.loc[m, "overlapping_gene_count"], errors="coerce")
    if v.empty:
        return 0
    mx = v.max()
    return int(mx) if pd.notna(mx) else 0


def _print_run_validation(
    *,
    datasets_used: list[str],
    used_mock: bool,
    pathway_scores: pd.DataFrame,
    top_genes: pd.DataFrame,
    mapping_loaded: bool,
    osd569_in_run: bool,
    osd569_pathway_overlap: bool,
    load_report: dict,
    narrative_path: Path,
    summary_json_path: Path,
    heatmap_path: Path,
    pathway_level_top_genes_export: pd.DataFrame | None = None,
) -> None:
    """End-of-run validation block (exploratory; does not fail the pipeline)."""
    print("\n=== Run validation ===")
    print(f"  Datasets analyzed: {', '.join(datasets_used)}")
    print(f"  Used mock data: {'yes' if used_mock else 'no'}")
    print(f"  Ensembl mapping loaded (usable pairs): {'yes' if mapping_loaded else 'no'}")
    print(
        f"  OSD-569 has pathway overlap (any overlapping genes > 0): "
        f"{'yes' if osd569_pathway_overlap else 'no'}"
    )

    if osd569_in_run:
        om = load_report.get("osd569_ensembl_mapping") or {}
        tu = int(om.get("total_unique_ensembl_ids", 0) or 0)
        mu = int(om.get("mapped_unique_ids", 0) or 0)
        z = _osd569_max_pathway_overlap(pathway_scores)
        if not mapping_loaded:
            print(
                "\n  OSD-569 loaded, but Ensembl-to-HGNC mapping was not loaded. "
                "OSD-569 pathway overlap is expected to be zero until "
                "data/processed/ensembl_to_symbol.tsv is added."
            )
        else:
            print(
                f"\n  OSD-569 mapping loaded. {mu} of {tu} unique Ensembl IDs mapped "
                f"to HGNC symbols. OSD-569 max pathway overlap after mapping: {z}."
            )

    if pathway_scores.empty:
        print("  All pathway score rows: 0")
        print("  Evidence-filtered pathway rows (overlap >= 3): 0")
        print("  Low-overlap pathway rows (overlap < 3): 0")
        print("  Top genes from pathway-level signals (overlap >= 3 parent): 0")
        print(f"  Narrative: {narrative_path}")
        print(f"  Summary JSON: {summary_json_path}")
        print(f"  Heatmap: {heatmap_path}")
        return

    ps = pathway_scores
    n_all = len(ps)
    n_ev = int((ps["overlapping_gene_count"].astype(float) >= 3).sum())
    n_low = int((ps["overlapping_gene_count"].astype(float) < 3).sum())
    if pathway_level_top_genes_export is not None:
        n_pl_genes = len(pathway_level_top_genes_export)
    else:
        pl_genes = pathway_level_top_genes_table(top_genes, load_report)
        n_pl_genes = len(pl_genes)

    print(f"  All pathway rows: {n_all}")
    print(f"  Evidence-filtered pathway rows (overlap >= 3): {n_ev}")
    print(f"  Low-overlap pathway rows (overlap < 3): {n_low}")
    print(f"  Gene rows (parent pathway overlap >= 3): {n_pl_genes}")

    if n_ev == 0:
        print(
            "\n  No pathway-level candidate rows found. "
            "Check gene identifiers, mapping, and gene-set overlap."
        )

    pl = ps[ps["overlapping_gene_count"].astype(float) >= 3].copy()
    if not pl.empty:
        best = pl.assign(
            _ps=pd.to_numeric(pl["pathway_rank_metric"], errors="coerce")
        ).sort_values("_ps", ascending=False, na_position="last").iloc[0]
        pm = best["pathway_rank_metric"]
        pm_fmt = f"{float(pm):.4f}" if pd.notna(pm) else "nan"
        print(
            "  Strongest pathway-level signal (overlap >= 3): "
            f"{best.get('pathway')} | {best.get('dataset_name')} | {best.get('cell_type')} | "
            f"pathway_rank_metric={pm_fmt}"
        )
    else:
        print("  Strongest pathway-level signal: (none — no overlap >= 3 rows)")

    gl = ps[ps["overlapping_gene_count"].astype(float).isin([1, 2])].copy()
    if not gl.empty:
        gr = gl.assign(
            _ps=pd.to_numeric(gl["pathway_rank_metric"], errors="coerce")
        ).sort_values("_ps", ascending=False, na_position="last").iloc[0]
        pm = gr["pathway_rank_metric"]
        pm_fmt = f"{float(pm):.4f}" if pd.notna(pm) else "nan"
        print(
            "  Strongest gene-level observation (overlap 1–2; not pathway evidence): "
            f"{gr.get('pathway')} | {gr.get('dataset_name')} | {gr.get('cell_type')} | "
            f"n={int(gr['overlapping_gene_count'])} | pathway_rank_metric={pm_fmt}"
        )
    else:
        print("  Strongest gene-level observation: (none)")

    if osd569_in_run and not osd569_pathway_overlap:
        if mapping_loaded:
            print(
                "\n  Mapping file loaded, but OSD-569 still has zero curated pathway overlap. "
                "Check whether the mapping symbols match the curated pathway gene symbols "
                "(see gene_sets.py and outputs/logs/data_debug_report.txt)."
            )
        else:
            print(
                "\n  OSD-569 has zero curated pathway overlap in this run. "
                "This is usually because the table uses Ensembl IDs and no mapping file was loaded."
            )

    print(f"  Narrative: {narrative_path}")
    print(f"  Summary JSON: {summary_json_path}")
    print(f"  Heatmap: {heatmap_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Exploratory DNA damage-response signature scoring from astronaut DE tables."
        )
    )
    parser.add_argument(
        "--download-from-notebook",
        action="store_true",
        help=(
            "Scan the Torchlight starter .ipynb for NASA OSDR URLs and download "
            "OSD-569 whole blood + OSD-570 PBMC gene-expression XLSX files into data/raw/."
        ),
    )
    parser.add_argument(
        "--notebook",
        type=Path,
        default=None,
        help=(
            "Path to Torchlight_Hackathon_2026.ipynb "
            "(default: data/raw/Torchlight_Hackathon_2026.ipynb)."
        ),
    )
    return parser.parse_args()


def standardized_column_schema() -> dict:
    """Describe columns after load_data standardization (for reporting)."""
    return {
        "gene": "HGNC symbol when available; otherwise Ensembl ID for OSD-569",
        "effect_size": "log2 fold change (dataset-specific contrast)",
        "adjusted_p_value": "adjusted p-value when available",
        "raw_p_value": "raw p-value when available",
        "comparison": "Pre-defined contrast label from each OSD study",
        "dataset_name": "OSD-569_whole_blood_RNA | OSD-570_PBMC_snRNA | MOCK*",
        "cell_type": "Immune subset for PBMC; 'unknown' for whole blood",
        "timepoint": "Reserved for future use ('unknown' in v1)",
    }


def main() -> int:
    args = parse_args()
    ensure_output_dirs()
    root = project_root()

    if args.download_from_notebook:
        from osdr_notebook_fetch import fetch_gene_expression_tables

        fetch_results = fetch_gene_expression_tables(
            root, notebook_path=args.notebook
        )
        print("\n=== OSDR fetch from Torchlight notebook ===")
        for key in sorted(fetch_results.keys()):
            print(f"  {key}: {fetch_results[key]}")
        print()

    try:
        combined_df, load_report, used_mock = load_all_tables(root)
    except Exception as exc:
        print("Fatal error while loading data:", exc)
        print(traceback.format_exc())
        return 1

    datasets_used = sorted(combined_df["dataset_name"].astype(str).unique().tolist())

    pathway_scores, top_genes = score_all_pathways(combined_df)
    pathway_scores = attach_interpretation_labels(pathway_scores)

    mapping_loaded = bool(load_report.get("mapping_loaded"))
    osd569_in_run = any("OSD-569" in str(d) for d in datasets_used)
    osd569_pathway_overlap = _osd569_has_any_pathway_overlap(pathway_scores)

    excluded_plot_datasets = datasets_without_pathway_identifier_overlap(combined_df)
    excluded_plot_datasets |= datasets_with_zero_pathway_overlap(pathway_scores)
    if excluded_plot_datasets:
        print("\n=== Plot exclusion (identifier mismatch or zero pathway overlap) ===")
        for ds in sorted(excluded_plot_datasets):
            print(
                f"  Dataset '{ds}' is omitted from pathway figures "
                "(no curated-symbol overlap and/or all pathway rows have overlap=0). "
                "See outputs/logs/data_debug_report.txt."
            )

    debug_report_path = write_data_debug_report(
        combined_df,
        load_report,
        root,
        pathway_scores=pathway_scores,
    )

    column_detection = {
        "standardized_schema": standardized_column_schema(),
        "heterogeneous_table_detection": (
            "Use column_detection.py on raw OSD exports when adding parsers; "
            "this run used standardized column names after study-specific loaders."
        ),
    }
    if load_report.get("osd570_raw_column_detection"):
        column_detection["osd570_raw_column_detection"] = load_report[
            "osd570_raw_column_detection"
        ]

    osd569_inclusion = build_osd569_score_inclusion_status(
        load_report,
        osd569_pathway_overlap=osd569_pathway_overlap,
        mapping_loaded=mapping_loaded,
    )

    pathway_level_export, osd569_integrity = prepare_pathway_level_top_genes_export(
        top_genes, load_report, root
    )
    if int(osd569_integrity.get("mismatched_rows", 0) or 0) > 0:
        print("\n=== WARNING: OSD-569 top-gene mapping integrity ===")
        print(
            f"  Rows failing Ensembl→gene validation (excluded from pathway-level export): "
            f"{osd569_integrity['mismatched_rows']}"
        )
        print("  See outputs/logs/osd569_mapping_integrity_check.txt")

    summary_json = build_initial_analysis_summary_json(
        datasets_used=datasets_used,
        column_detection=column_detection,
        pathway_scores=pathway_scores,
        top_genes=top_genes,
        load_report=load_report,
        used_mock=used_mock,
        mapping_loaded=mapping_loaded,
        osd569_has_pathway_overlap=osd569_pathway_overlap,
        osd569_score_inclusion_status=osd569_inclusion,
        pathway_level_top_genes_df=pathway_level_export,
        osd569_mapping_integrity_status=osd569_integrity,
    )

    save_summary_outputs(
        pathway_scores,
        top_genes,
        summary_json,
        root,
        load_report=load_report,
        pathway_level_top_genes_df=pathway_level_export,
    )
    summary_json_path = root / "outputs" / "logs" / "initial_analysis_summary.json"
    heatmap_path = root / "outputs" / "figures" / "evidence_filtered_pathway_heatmap.png"
    narrative_path = build_analysis_narrative_md(
        pathway_scores=pathway_scores,
        load_report=load_report,
        datasets_used=datasets_used,
        used_mock=used_mock,
        osd569_has_pathway_overlap=osd569_pathway_overlap,
        mapping_loaded=mapping_loaded,
        top_genes=top_genes,
        osd569_numeric_status=load_report.get("osd569_numeric_status"),
        pathway_level_top_genes_df=pathway_level_export,
        osd569_mapping_integrity=osd569_integrity,
    )

    hm_extra: set[str] = set()
    hm_footer: str | None = None
    ons = load_report.get("osd569_numeric_status") or {}
    if (
        osd569_in_run
        and ons
        and not bool(ons.get("osd569_scoreable"))
        and osd569_pathway_overlap
    ):
        hm_extra.add("OSD-569_whole_blood_RNA")
        hm_footer = (
            "OSD-569 overlap found, but numeric effect sizes were not parsed."
        )

    run_all_plots(
        pathway_scores,
        top_genes,
        excluded_plot_datasets,
        root,
        heatmap_extra_exclusions=hm_extra,
        heatmap_footer_warning=hm_footer,
    )

    rep = ColumnDetectionReport(
        gene_column="gene",
        effect_column="effect_size",
        adjusted_p_column="adjusted_p_value",
        raw_p_column="raw_p_value",
        meta_columns={
            "comparison": "comparison",
            "cell_type": "cell_type",
            "dataset_name": "dataset_name",
        },
        notes=[
            "Post-load standardized schema; see references/PROVENANCE.md for raw OSD layouts.",
        ],
    )

    print("\n=== Column detection (standardized table) ===")
    print(rep.to_dict())

    print("\n=== Datasets run ===")
    print(", ".join(datasets_used))
    if used_mock:
        print("(Mock data — place official XLSX files in data/raw/ for real OSD runs.)")

    _print_run_validation(
        datasets_used=datasets_used,
        used_mock=used_mock,
        pathway_scores=pathway_scores,
        top_genes=top_genes,
        mapping_loaded=mapping_loaded,
        osd569_in_run=osd569_in_run,
        osd569_pathway_overlap=osd569_pathway_overlap,
        load_report=load_report,
        narrative_path=narrative_path,
        summary_json_path=summary_json_path,
        heatmap_path=heatmap_path,
        pathway_level_top_genes_export=pathway_level_export,
    )

    inspect_wb = root / "outputs" / "logs" / "osd569_workbook_inspection.txt"
    inspect_hint = ""
    if not inspect_wb.is_file():
        inspect_hint = (
            "\nFor workbook structure details, run: python scripts/inspect_osd569_workbook.py\n"
        )

    print(
        f"\nDebug report: {debug_report_path}\n"
        "Outputs written under outputs/tables, outputs/figures, outputs/logs."
        f"{inspect_hint}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
