"""
Compact summaries for CSV exports and initial_analysis_summary.json (AutoBioScout hook).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from gene_id_mapping import load_ensembl_symbol_mapping, strip_ensembl_version
from utils import project_root, write_csv, write_json


def pathway_numeric_evidence_mask(ps: pd.DataFrame) -> pd.Series:
    """True where pathway_rank_metric and average_signed_effect are both numeric."""
    prm = pd.to_numeric(ps["pathway_rank_metric"], errors="coerce")
    ase = pd.to_numeric(ps["average_signed_effect"], errors="coerce")
    return prm.notna() & ase.notna()


def split_overlap_evidence_tables(ps: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    overlap ≥3 with numeric scores vs overlap ≥3 but missing pathway numeric metrics.
    """
    ov = ps["overlapping_gene_count"].astype(float) >= 3
    nm = pathway_numeric_evidence_mask(ps)
    return ps[ov & nm].copy(), ps[ov & ~nm].copy()


def pathway_scores_for_summary_json(
    pathway_scores: pd.DataFrame,
    load_report: dict[str, Any] | None,
) -> pd.DataFrame:
    """Exclude OSD-569 rows unless scoreable and numerically evidenced."""
    lr = load_report or {}
    num = lr.get("osd569_numeric_status") or {}
    sc = bool(num.get("osd569_scoreable"))
    is569 = pathway_scores["dataset_name"].astype(str).str.contains("OSD-569", na=False)
    nm = pathway_numeric_evidence_mask(pathway_scores)
    keep = (~is569 & nm) | (is569 & sc & nm)
    return pathway_scores.loc[keep].copy()


def build_osd569_score_inclusion_status(
    load_report: dict[str, Any] | None,
    *,
    osd569_pathway_overlap: bool,
    mapping_loaded: bool,
) -> dict[str, Any]:
    """Summarize whether OSD-569 rows belong in score tables / heatmap."""
    num = (load_report or {}).get("osd569_numeric_status") or {}
    sc = bool(num.get("osd569_scoreable"))
    has_eff = int(num.get("effect_size_non_null_count") or 0) > 0
    reason = ""
    if not mapping_loaded:
        reason = "Ensembl mapping was not loaded or had zero usable pairs."
    elif not has_eff:
        reason = "No numeric effect sizes were parsed from the OSD-569 workbook."
    elif not sc:
        reason = (
            "OSD-569 did not meet scoreability criteria "
            "(requires mapping loaded, non-null effect sizes, and mapped curated genes with numeric effects)."
        )
    return {
        "mapping_loaded": mapping_loaded,
        "has_pathway_overlap": osd569_pathway_overlap,
        "has_numeric_effect_sizes": has_eff,
        "included_in_score_tables": sc,
        "included_in_heatmap": sc,
        "reason_if_excluded": reason if not sc else "",
    }


def filter_top_genes_for_score_summaries(
    top_genes: pd.DataFrame,
    load_report: dict[str, Any] | None,
) -> pd.DataFrame:
    """
    Drop rows with NaN gene_signal/effect_size; drop all OSD-569 rows when not scoreable.
    Mapped OSD-569 genes are excluded from top-gene summaries when numeric parsing failed.
    """
    if top_genes.empty:
        return top_genes
    tg = top_genes.copy()
    lr = load_report or {}
    num = lr.get("osd569_numeric_status") or {}
    sc = bool(num.get("osd569_scoreable"))
    is569 = tg["dataset_name"].astype(str).str.contains("OSD-569", na=False)
    bad = pd.to_numeric(tg["gene_signal"], errors="coerce").isna() | pd.to_numeric(
        tg["effect_size"], errors="coerce"
    ).isna()
    drop = bad | (is569 & (~pd.Series([sc] * len(tg), index=tg.index)))
    return tg.loc[~drop].copy()


def build_analysis_summary_csv(pathway_scores: pd.DataFrame, load_report: dict[str, Any] | None = None) -> pd.DataFrame:
    """Highest exploratory ranking metrics — pathway-level rows (overlap ≥3) only."""
    ps = pathway_scores_for_summary_json(pathway_scores, load_report)
    pl = ps[ps["overlapping_gene_count"].astype(float) >= 3].copy()
    if pl.empty:
        return pd.DataFrame(
            columns=[
                "rank",
                "dataset_name",
                "comparison",
                "cell_type",
                "pathway",
                "pathway_rank_metric",
                "interpretation_label",
            ]
        )

    sortable = pl.copy()
    sortable["_score"] = pd.to_numeric(
        sortable["pathway_rank_metric"], errors="coerce"
    )
    sortable = sortable.sort_values("_score", ascending=False)
    sortable = sortable.drop(columns=["_score"])
    sortable.insert(0, "rank", range(1, len(sortable) + 1))
    return sortable.head(25)


def _top_contributing_genes_row(
    row: pd.Series,
    top_genes: pd.DataFrame,
    max_genes: int = 10,
) -> str:
    sub = top_genes[
        (top_genes["dataset_name"].astype(str) == str(row["dataset_name"]))
        & (top_genes["comparison"].astype(str) == str(row["comparison"]))
        & (top_genes["cell_type"].astype(str) == str(row["cell_type"]))
        & (top_genes["pathway"].astype(str) == str(row["pathway"]))
    ].copy()
    gs_ok = pd.to_numeric(sub["gene_signal"], errors="coerce").notna()
    ef_ok = pd.to_numeric(sub["effect_size"], errors="coerce").notna()
    sub = sub.loc[gs_ok & ef_ok].sort_values("gene_signal", ascending=False)
    sub = sub.head(max_genes)
    if sub.empty:
        return ""
    return "; ".join(sub["gene"].astype(str).tolist())


PATHWAY_LEVEL_GENE_COLUMNS = [
    "gene",
    "original_gene_id",
    "ensembl_id",
    "pathway",
    "effect_size",
    "raw_p_value",
    "adjusted_p_value",
    "gene_signal",
    "direction",
    "dataset_name",
    "comparison",
    "cell_type",
    "evidence_tier_of_parent_pathway",
    "parent_pathway_overlap_count",
]


# Regression spot-checks for OSD-569 Ensembl → symbol consistency (export table).
OSD569_CANONICAL_ENSEMBL_GENE: dict[str, str] = {
    "ENSG00000012048": "BRCA1",
    "ENSG00000079246": "XRCC5",
    "ENSG00000012061": "ERCC1",
    "ENSG00000141510": "TP53",
    "ENSG00000149311": "ATM",
}


def filter_main_candidate_pathway_rows(
    pathway_scores: pd.DataFrame,
    load_report: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """
    Rows for main_candidate_signals.csv and narrative "main pathway-level candidates".

    overlap >= 3 AND numeric pathway scores AND (pathway_rank_metric >= 0.5 OR significant_gene_count >= 1).

    When OSD-569 is not scoreable (numeric DE / mapping criteria), **no** OSD-569 rows are returned
    (acceptance test SUCCESS B).
    """
    ps = pathway_scores[pathway_scores["overlapping_gene_count"].astype(float) >= 3].copy()
    if ps.empty:
        return ps
    ps = ps[pathway_numeric_evidence_mask(ps)].copy()
    if ps.empty:
        return ps
    lr = load_report or {}
    num = lr.get("osd569_numeric_status") or {}
    if not bool(num.get("osd569_scoreable")):
        is569 = ps["dataset_name"].astype(str).str.contains("OSD-569", na=False)
        ps = ps.loc[~is569].copy()
    if ps.empty:
        return ps
    prm = pd.to_numeric(ps["pathway_rank_metric"], errors="coerce")
    sig = pd.to_numeric(ps["significant_gene_count"], errors="coerce").fillna(0)
    mask = (prm >= 0.5) | (sig >= 1)
    return ps.loc[mask].copy()


def pathway_level_top_genes_table(
    top_genes: pd.DataFrame,
    load_report: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """
    Gene rows whose parent pathway has overlap >= 3 (pathway-level evidence only).
    """
    if top_genes.empty:
        return pd.DataFrame(columns=PATHWAY_LEVEL_GENE_COLUMNS)
    tg = filter_top_genes_for_score_summaries(top_genes, load_report)
    if tg.empty:
        return pd.DataFrame(columns=PATHWAY_LEVEL_GENE_COLUMNS)
    if "parent_pathway_overlap_count" not in tg.columns:
        return pd.DataFrame(columns=PATHWAY_LEVEL_GENE_COLUMNS)
    poc = pd.to_numeric(tg["parent_pathway_overlap_count"], errors="coerce").fillna(0)
    out = tg.loc[poc >= 3].copy()
    for c in PATHWAY_LEVEL_GENE_COLUMNS:
        if c not in out.columns:
            out[c] = pd.NA
    out = out[PATHWAY_LEVEL_GENE_COLUMNS]
    out = out.loc[
        pd.to_numeric(out["gene_signal"], errors="coerce").notna()
        & pd.to_numeric(out["effect_size"], errors="coerce").notna()
    ].copy()
    return out.sort_values(
        ["dataset_name", "comparison", "cell_type", "pathway", "gene_signal"],
        ascending=[True, True, True, True, False],
    )


def _osd569_canonical_regression_checks(pl_export: pd.DataFrame) -> list[dict[str, Any]]:
    """Spot-check known Ensembl IDs vs gene symbols on the final export table."""
    checks: list[dict[str, Any]] = []
    sub = pl_export[pl_export["dataset_name"].astype(str) == "OSD-569_whole_blood_RNA"].copy()
    for ens_base, sym_exp in OSD569_CANONICAL_ENSEMBL_GENE.items():
        sel = sub[
            sub["original_gene_id"].map(lambda x: strip_ensembl_version(str(x)) == ens_base)
        ]
        present = not sel.empty
        rec: dict[str, Any] = {
            "ensembl_id": ens_base,
            "expected_symbol": sym_exp,
            "present_in_export": present,
        }
        if present:
            ok = True
            for _, rr in sel.iterrows():
                if str(rr.get("gene", "")).strip().upper() != sym_exp.upper():
                    ok = False
                    break
            rec["gene_matches_expected"] = ok
            rec["example_observed_gene"] = str(sel.iloc[0].get("gene", ""))
        else:
            rec["gene_matches_expected"] = None
            rec["example_observed_gene"] = None
        checks.append(rec)
    return checks


def prepare_pathway_level_top_genes_export(
    top_genes: pd.DataFrame,
    load_report: dict[str, Any] | None,
    root: Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Build pathway-level top-gene table, drop OSD-569 rows that fail Ensembl→gene validation,
    write osd569_mapping_integrity_check.txt (including canonical regression checks).
    """
    root = root or project_root()
    logs = root / "outputs" / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    detail_path = logs / "osd569_mapping_integrity_check.txt"
    lr = load_report or {}

    status: dict[str, Any] = {
        "checked_rows": 0,
        "mismatched_rows": 0,
        "passed": True,
        "osd569_rows_removed": 0,
        "regression_checks": [],
        "canonical_regression_passed": True,
    }

    pl = pathway_level_top_genes_table(top_genes, lr)

    def _stamp_canonical() -> None:
        rc = status.get("regression_checks") or []
        status["canonical_regression_passed"] = all(
            not (r.get("present_in_export") and r.get("gene_matches_expected") is False)
            for r in rc
        )

    def _write_skip(msg: str) -> tuple[pd.DataFrame, dict[str, Any]]:
        detail_path.write_text(msg + "\n", encoding="utf-8")
        status["regression_checks"] = _osd569_canonical_regression_checks(pl)
        _stamp_canonical()
        return pl, status

    if pl.empty or "original_gene_id" not in pl.columns:
        return _write_skip("OSD-569 mapping integrity: skipped (empty pathway-level top genes or no original_gene_id).")

    sub569 = pl[pl["dataset_name"].astype(str) == "OSD-569_whole_blood_RNA"].copy()
    status["checked_rows"] = int(len(sub569))

    if sub569.empty:
        detail_lines = [
            "OSD-569 mapping integrity check",
            "(pathway-level top genes export)",
            "",
            "checked_rows (OSD-569 in export): 0",
            "mismatched_rows: 0",
            "passed: True",
            "osd569_rows_removed: 0",
            "",
        ]
        pl_clean = pl.copy()
        status["regression_checks"] = _osd569_canonical_regression_checks(pl_clean)
        detail_lines.append("### Canonical Ensembl regression (spot-checks)")
        for r in status["regression_checks"]:
            detail_lines.append(f"  {r['ensembl_id']} → expected {r['expected_symbol']!r} | "
                                f"present={r['present_in_export']} | match={r.get('gene_matches_expected')}")
        _stamp_canonical()
        detail_path.write_text("\n".join(detail_lines), encoding="utf-8")
        return pl_clean, status

    if not lr.get("mapping_loaded"):
        status["regression_checks"] = _osd569_canonical_regression_checks(pl)
        _stamp_canonical()
        detail_path.write_text(
            "OSD-569 mapping integrity: skipped (Ensembl→symbol mapping not loaded).\n",
            encoding="utf-8",
        )
        return pl, status

    mp = lr.get("mapping_path")
    path = Path(str(mp)) if mp else Path()
    if not path.is_file():
        status["regression_checks"] = _osd569_canonical_regression_checks(pl)
        _stamp_canonical()
        detail_path.write_text(
            f"OSD-569 mapping integrity: skipped (mapping file missing: {path}).\n",
            encoding="utf-8",
        )
        return pl, status

    mapping, _meta = load_ensembl_symbol_mapping(path)
    drop_ix: list[Any] = []
    mismatches: list[dict[str, Any]] = []

    for ix, row in sub569.iterrows():
        oid = row.get("original_gene_id")
        cleaned = strip_ensembl_version(str(oid))
        exp = mapping.get(cleaned)
        if exp is None or not str(exp).strip():
            continue
        obs = str(row.get("gene", "")).strip()
        if str(exp).strip().upper() != obs.upper():
            drop_ix.append(ix)
            mismatches.append(
                {
                    "original_gene_id": oid,
                    "cleaned_ensembl_id": cleaned,
                    "expected_symbol_from_mapping": exp,
                    "observed_gene_symbol": obs,
                    "pathway": row.get("pathway"),
                    "effect_size": row.get("effect_size"),
                    "adjusted_p_value": row.get("adjusted_p_value"),
                    "gene_signal": row.get("gene_signal"),
                }
            )

    status["mismatched_rows"] = len(mismatches)
    status["passed"] = len(mismatches) == 0

    pl_clean = pl.drop(index=drop_ix) if drop_ix else pl.copy()
    status["osd569_rows_removed"] = len(drop_ix)

    status["regression_checks"] = _osd569_canonical_regression_checks(pl_clean)

    detail_lines: list[str] = [
        "OSD-569 mapping integrity check",
        "(pathway-level top genes — after row-identity fix and optional row drops)",
        "",
        f"checked_rows (OSD-569 rows before drop): {status['checked_rows']}",
        f"mismatched_rows (failed Ensembl→gene vs mapping file): {status['mismatched_rows']}",
        f"passed (no mismatches among verifiable rows): {status['passed']}",
        f"osd569_rows_removed from export: {status['osd569_rows_removed']}",
        "",
    ]

    if mismatches:
        detail_lines.append("Mismatch detail (rows excluded from export):")
        detail_lines.append(pd.DataFrame(mismatches).to_string(index=False))
        detail_lines.append("")
        dbg_append = logs / "data_debug_report.txt"
        if dbg_append.is_file():
            with dbg_append.open("a", encoding="utf-8") as f:
                f.write("\n### OSD-569 mapping integrity (top genes)\n\n")
                f.write(
                    "WARNING: One or more OSD-569 pathway-level top-gene rows disagreed with "
                    "ensembl_to_symbol.tsv; those rows were excluded from the export.\n"
                )
                f.write(
                    f"See outputs/logs/osd569_mapping_integrity_check.txt ({len(mismatches)} row(s)).\n"
                )
    elif status["checked_rows"] > 0:
        detail_lines.append(
            "OSD-569 mapping integrity check passed: all verifiable top-gene rows match "
            "original_gene_id → gene symbol mapping."
        )
        detail_lines.append("")

    detail_lines.append("### Canonical Ensembl regression (spot-checks on final export)")
    for r in status["regression_checks"]:
        pres = r["present_in_export"]
        match = r.get("gene_matches_expected")
        ex = r.get("example_observed_gene")
        detail_lines.append(
            f"  {r['ensembl_id']} → expected {r['expected_symbol']!r} | "
            f"present={pres} | match={match} | observed={ex!r}"
        )

    _stamp_canonical()

    detail_path.write_text("\n".join(detail_lines), encoding="utf-8")
    return pl_clean, status


def build_initial_analysis_summary_json(
    *,
    datasets_used: list[str],
    column_detection: dict[str, Any],
    pathway_scores: pd.DataFrame,
    top_genes: pd.DataFrame,
    load_report: dict[str, Any],
    used_mock: bool,
    mapping_loaded: bool,
    osd569_has_pathway_overlap: bool,
    osd569_score_inclusion_status: dict[str, Any] | None = None,
    pathway_level_top_genes_df: pd.DataFrame | None = None,
    osd569_mapping_integrity_status: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Structured summary — separates pathway-level vs gene-level evidence."""

    limitations = [
        "This workflow summarizes exploratory DNA damage-response–associated "
        "gene expression signatures; it does not detect DNA breaks or lesions.",
        "The Inspiration4 crew sample size is small (n=4); statistical "
        "generalization is limited.",
        "Transcriptional shifts may reflect inflammation, stress, immune "
        "reorganization, or cell-composition changes—not only DNA repair biology.",
        "Results are exploratory and not diagnostic or clinical.",
    ]

    ps = pathway_scores_for_summary_json(pathway_scores, load_report)
    _bad_tiers = {"gene-level signal only", "no overlapping genes"}
    pathway_level = ps[
        (ps["overlapping_gene_count"].astype(float) >= 3)
        & (~ps["evidence_tier"].astype(str).str.strip().isin(_bad_tiers))
    ].copy()
    pathway_level = pathway_level[pathway_numeric_evidence_mask(pathway_level)].copy()
    gene_level_only = ps[ps["overlapping_gene_count"].astype(float).isin([1, 2])].copy()

    strongest_pathway_level_signal = None
    if not pathway_level.empty:
        row = pathway_level.assign(
            _ps=pd.to_numeric(pathway_level["pathway_rank_metric"], errors="coerce")
        ).sort_values("_ps", ascending=False, na_position="last").iloc[0]
        strongest_pathway_level_signal = {
            "pathway": row.get("pathway"),
            "dataset_name": row.get("dataset_name"),
            "comparison": row.get("comparison"),
            "cell_type": row.get("cell_type"),
            "pathway_rank_metric": row.get("pathway_rank_metric"),
            "interpretation_label": row.get("interpretation_label"),
            "evidence_tier": row.get("evidence_tier"),
            "overlapping_gene_count": int(row.get("overlapping_gene_count", 0) or 0),
        }

    strongest_gene_level_signal = None
    if not gene_level_only.empty:
        row = gene_level_only.assign(
            _ps=pd.to_numeric(gene_level_only["pathway_rank_metric"], errors="coerce")
        ).sort_values("_ps", ascending=False, na_position="last").iloc[0]
        strongest_gene_level_signal = {
            "pathway": row.get("pathway"),
            "dataset_name": row.get("dataset_name"),
            "comparison": row.get("comparison"),
            "cell_type": row.get("cell_type"),
            "pathway_rank_metric": row.get("pathway_rank_metric"),
            "evidence_tier": row.get("evidence_tier"),
            "overlapping_gene_count": int(row.get("overlapping_gene_count", 0) or 0),
            "note": (
                "Gene-level observation only — not sufficient for pathway-level interpretation."
            ),
        }

    highest_candidate_pathway_signals: list[str] = []
    main_candidate_rows = filter_main_candidate_pathway_rows(ps, load_report)
    if not main_candidate_rows.empty:
        hi = main_candidate_rows.sort_values(
            "pathway_rank_metric", ascending=False, na_position="last"
        )
        for _, r in hi.head(15).iterrows():
            prm = pd.to_numeric(r.get("pathway_rank_metric"), errors="coerce")
            prm_fmt = f"{float(prm):.3f}" if pd.notna(prm) else "nan"
            highest_candidate_pathway_signals.append(
                f"{r.get('dataset_name')} | {r.get('comparison')} | {r.get('cell_type')} | "
                f"{r.get('pathway')}: candidate pathway-level exploratory signal "
                f"(rank_metric={prm_fmt}, "
                f"n={int(r['overlapping_gene_count'])})"
            )

    highest_gene_level_observations: list[str] = []
    if not gene_level_only.empty:
        sg = gene_level_only.sort_values(
            "pathway_rank_metric", ascending=False, na_position="last"
        )
        for _, r in sg.head(15).iterrows():
            highest_gene_level_observations.append(
                f"{r.get('dataset_name')} | {r.get('cell_type')} | {r.get('pathway')}: "
                f"gene-level observation only (n={int(r['overlapping_gene_count'])}) — "
                "not pathway-level evidence."
            )

    tg_sum = filter_top_genes_for_score_summaries(top_genes, load_report)

    strongest_genes: list[str] = []
    if not tg_sum.empty:
        strongest_genes = (
            tg_sum.sort_values("gene_signal", ascending=False)
            .head(15)["gene"]
            .astype(str)
            .tolist()
        )

    pl_genes = (
        pathway_level_top_genes_df
        if pathway_level_top_genes_df is not None
        else pathway_level_top_genes_table(top_genes, load_report)
    )
    pathway_level_gene_hits: list[str] = []
    if not pl_genes.empty:
        for _, r in pl_genes.head(30).iterrows():
            gs = r.get("gene_signal")
            if pd.isna(pd.to_numeric(gs, errors="coerce")):
                continue
            pathway_level_gene_hits.append(
                f"{r.get('gene')} | {r.get('pathway')} | {r.get('dataset_name')} | "
                f"{r.get('cell_type')} | gene_signal={float(gs):.4f}"
            )

    osd569_pathway_interpretation_note: str | None = None
    if "OSD-569_whole_blood_RNA" in datasets_used and not osd569_has_pathway_overlap:
        if not mapping_loaded:
            osd569_pathway_interpretation_note = (
                "OSD-569 still has no overlap with curated pathway symbols in this run because "
                "the table uses Ensembl IDs and no successful Ensembl-to-HGNC mapping was loaded. "
                "Until mapping is added, this pipeline cannot determine whether OSD-569 contains "
                "DNA-damage-response pathway genes."
            )
        else:
            osd569_pathway_interpretation_note = (
                "Mapping file loaded, but OSD-569 still has zero curated pathway overlap in this run. "
                "Check whether mapped symbols align with curated pathway gene symbols (gene_sets.py), "
                "mapping coverage, and whether differential-expression genes intersect those lists."
            )

    overlaps: dict[str, int] = {}
    for _, r in pathway_scores.iterrows():
        key = (
            f'{r.get("dataset_name", "")}::{r.get("comparison", "")}::'
            f'{r.get("cell_type", "")}::{r.get("pathway", "")}'
        )
        overlaps[key] = int(r.get("overlapping_gene_count", 0) or 0)

    osd569_numeric_block = (load_report or {}).get("osd569_numeric_status")

    osd569_top_gene_exclusion_note = None
    if (
        osd569_numeric_block
        and not osd569_numeric_block.get("osd569_scoreable")
        and "OSD-569_whole_blood_RNA" in datasets_used
    ):
        osd569_top_gene_exclusion_note = (
            "Mapped OSD-569 genes were excluded from top-gene summaries because numeric effect "
            "sizes were not parsed."
        )

    summary: dict[str, Any] = {
        "project": "Spaceflight DNA Damage Response Classifier",
        "datasets_analyzed": datasets_used,
        "used_mock_data": used_mock,
        "ensembl_mapping_loaded": mapping_loaded,
        "osd569_has_any_pathway_overlap": osd569_has_pathway_overlap,
        "osd569_numeric_status": osd569_numeric_block,
        "osd569_score_inclusion_status": osd569_score_inclusion_status,
        "columns_detected_or_standardized": column_detection,
        "strongest_pathway_level_signal": strongest_pathway_level_signal,
        "strongest_gene_level_signal": strongest_gene_level_signal,
        "highest_candidate_pathway_signals": highest_candidate_pathway_signals,
        "highest_gene_level_observations": highest_gene_level_observations,
        "strongest_genes_ranked_by_gene_signal": strongest_genes,
        "pathway_level_gene_hits": pathway_level_gene_hits,
        "overlapping_genes_per_pathway_row": overlaps,
        "limitations": limitations,
        "ranking_metric_note": (
            "pathway_rank_metric combines effect size and significance weight for exploratory "
            "ranking — not clinical severity or DNA lesion burden."
        ),
        "suggested_follow_up_analyses": [
            "Expand OSD-569 from curated-only Ensembl→HGNC mapping to a genome-wide annotation map if broader pathway libraries are added.",
            "Integrate OSD-570 outputs across cell types with composition-aware models.",
            "Compare concordance between whole blood and PBMC immune subsets.",
            "Validate hits with orthogonal assays (not inferred from RNA alone).",
            "Prefer clearer contrasts when available: flight perturbation (e.g. R+1 vs pre-flight) vs recovery shift (R+45 vs R+1).",
        ],
        "load_report": load_report,
    }

    if osd569_pathway_interpretation_note:
        summary["osd569_pathway_interpretation_note"] = osd569_pathway_interpretation_note

    if osd569_top_gene_exclusion_note:
        summary["osd569_top_gene_exclusion_note"] = osd569_top_gene_exclusion_note

    if osd569_mapping_integrity_status is not None:
        summary["osd569_mapping_integrity_status"] = osd569_mapping_integrity_status
        mi = osd569_mapping_integrity_status
        if int(mi.get("mismatched_rows", 0) or 0) > 0:
            summary["osd569_top_gene_mapping_integrity_note"] = (
                "OSD-569 top-gene rows were excluded because original_gene_id to gene-symbol "
                "integrity validation failed (see outputs/logs/osd569_mapping_integrity_check.txt)."
            )

    return summary


PATHWAY_SCORE_COLUMNS = [
    "dataset_name",
    "comparison",
    "cell_type",
    "pathway",
    "pathway_rank_metric",
    "evidence_tier",
    "average_signed_effect",
    "average_absolute_effect",
    "overlapping_gene_count",
    "upregulated_gene_count",
    "downregulated_gene_count",
    "significant_gene_count",
    "direction",
    "interpretation_label",
]

MAIN_CANDIDATE_COLUMNS = PATHWAY_SCORE_COLUMNS + ["top_contributing_genes"]


def save_summary_outputs(
    pathway_scores: pd.DataFrame,
    top_genes: pd.DataFrame,
    summary_payload: dict[str, Any],
    root: Path | None = None,
    load_report: dict[str, Any] | None = None,
    pathway_level_top_genes_df: pd.DataFrame | None = None,
) -> None:
    root = root or project_root()
    lr = load_report or {}
    ps = pathway_scores.copy()
    for c in PATHWAY_SCORE_COLUMNS:
        if c not in ps.columns:
            ps[c] = pd.NA
    ps = ps[PATHWAY_SCORE_COLUMNS]

    write_csv(ps, root / "outputs" / "tables" / "pathway_scores_all_rows.csv")

    ps_ev, ps_unscored = split_overlap_evidence_tables(ps)
    write_csv(ps_ev, root / "outputs" / "tables" / "pathway_scores_evidence_filtered.csv")
    write_csv(
        ps_unscored,
        root / "outputs" / "tables" / "pathway_scores_overlap_but_unscored.csv",
    )
    write_csv(ps_ev, root / "outputs" / "tables" / "pathway_scores.csv")

    low = ps[ps["overlapping_gene_count"].astype(float) < 3].copy()
    write_csv(low, root / "outputs" / "tables" / "low_overlap_signals.csv")

    # Main candidates: overlap ≥3 plus rank/significance gate (weak rows stay in ps_ev only).
    cand = filter_main_candidate_pathway_rows(ps, lr)
    tg_filt = filter_top_genes_for_score_summaries(top_genes, lr)
    cand["top_contributing_genes"] = cand.apply(
        lambda row: _top_contributing_genes_row(row, tg_filt),
        axis=1,
    )
    for c in MAIN_CANDIDATE_COLUMNS:
        if c not in cand.columns:
            cand[c] = pd.NA
    cand = cand[MAIN_CANDIDATE_COLUMNS]
    write_csv(cand, root / "outputs" / "tables" / "main_candidate_signals.csv")

    tg = top_genes.copy()
    write_csv(tg, root / "outputs" / "tables" / "top_damage_response_genes.csv")

    pl_tg = (
        pathway_level_top_genes_df
        if pathway_level_top_genes_df is not None
        else pathway_level_top_genes_table(top_genes, lr)
    )
    write_csv(
        pl_tg,
        root / "outputs" / "tables" / "top_genes_from_pathway_level_signals.csv",
    )

    write_csv(
        build_analysis_summary_csv(pathway_scores, lr),
        root / "outputs" / "tables" / "analysis_summary.csv",
    )
    write_json(summary_payload, root / "outputs" / "logs" / "initial_analysis_summary.json")
