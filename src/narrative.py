"""
Auto-generated markdown narrative for reviewers (exploratory, non-clinical).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from summarize_results import filter_main_candidate_pathway_rows, pathway_level_top_genes_table
from utils import project_root


def build_analysis_narrative_md(
    *,
    pathway_scores: pd.DataFrame,
    load_report: dict[str, Any],
    datasets_used: list[str],
    used_mock: bool,
    osd569_has_pathway_overlap: bool,
    mapping_loaded: bool,
    top_genes: pd.DataFrame,
    osd569_numeric_status: dict[str, Any] | None = None,
    pathway_level_top_genes_df: pd.DataFrame | None = None,
    osd569_mapping_integrity: dict[str, Any] | None = None,
) -> Path:
    """Write outputs/logs/analysis_narrative.md."""
    root = project_root()
    path = root / "outputs" / "logs" / "analysis_narrative.md"

    pl = filter_main_candidate_pathway_rows(pathway_scores, load_report)
    gl = pathway_scores[
        pathway_scores["overlapping_gene_count"].astype(float).isin([1, 2])
    ].copy()

    lines: list[str] = [
        "# Spaceflight DNA Damage Response Classifier — Initial Analysis Narrative",
        "",
        "**Opening warning:** Exploratory transcriptional signatures only. Does not detect DNA breaks, lesions, or clinical injury.",
        "",
        "## What data was analyzed",
        "",
        f"Datasets in this run: **{', '.join(datasets_used)}**.",
        "Targets are **OSD-569** whole-blood RNA and **OSD-570** PBMC snRNA when official XLSX files are placed in `data/raw/`.",
        "",
        f"Mock demonstration mode: **{'yes' if used_mock else 'no'}**.",
        "",
        "## What worked",
        "",
        "- **OSD-570 PBMC snRNA** uses HGNC-like gene symbols in the processed table, so overlap with curated pathway lists is computable without Ensembl mapping.",
        "",
    ]

    score569 = osd569_numeric_status.get("osd569_scoreable") if osd569_numeric_status else None

    if (
        "OSD-569_whole_blood_RNA" in datasets_used
        and mapping_loaded
        and score569 is True
    ):
        lines.extend(
            [
                "- **OSD-569:** Gene-symbol mapping and numeric differential-expression parsing succeeded, so OSD-569 "
                "contributes pathway-level whole-blood signals where overlap thresholds are met.",
                "",
                "OSD-569 now contributes scoreable whole-blood pathway rows with stronger gene-set coverage, but the "
                "observed mean signed effects are small and no OSD-569 pathway currently passes the main-candidate "
                "threshold. This summary is exploratory and does not describe clinical severity.",
                "",
            ]
        )

    needs_fixing: list[str] = []

    if "OSD-569_whole_blood_RNA" in datasets_used:
        if not mapping_loaded:
            needs_fixing.extend(
                [
                    "- **OSD-569** requires **`data/processed/ensembl_to_symbol.tsv`** (see `docs/ENSEMBL_MAPPING.md`). "
                    "Version suffixes on Ensembl IDs are stripped automatically.",
                    "",
                ]
            )
        if (
            mapping_loaded
            and osd569_has_pathway_overlap
            and score569 is False
        ):
            needs_fixing.append(
                "- **OSD-569:** Gene-symbol mapping succeeded and curated pathway overlap was found, but numeric "
                "differential-expression values were not parsed. OSD-569 is therefore excluded from score-based "
                "pathway summaries until the numeric loader is fixed (see `outputs/logs/osd569_loader_debug.txt`)."
            )
            needs_fixing.append("")
        if not osd569_has_pathway_overlap:
            if not mapping_loaded:
                needs_fixing.append(
                    "- **OSD-569 still has no overlap with curated pathway symbols in this run because the table uses "
                    "Ensembl IDs and no successful Ensembl-to-HGNC mapping was loaded. Until mapping is added, this "
                    "pipeline cannot determine whether OSD-569 contains DNA-damage-response pathway genes.**"
                )
            else:
                needs_fixing.append(
                    "- **OSD-569** shows **no pathway-level overlap** in this run after applying the mapping table "
                    "(see `data_debug_report.txt` for coverage — mapping completeness and DE gene content may still "
                    "limit overlap)."
                )
            needs_fixing.append("")

    lines.extend(
        [
            "## What needs fixing",
            "",
        ]
    )

    if needs_fixing:
        lines.extend(needs_fixing)
    elif (
        "OSD-569_whole_blood_RNA" in datasets_used
        and mapping_loaded
        and score569 is True
    ):
        lines.extend(
            [
                "No current loader blocker. OSD-569 is now scoreable, but its whole-blood effects are small and do "
                "not pass the main-candidate threshold.",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "_No blocking loader issues identified for the datasets in this run._",
                "",
            ]
        )

    lines.extend(
        [
            "## Main pathway-level candidate signals",
            "",
            "*Rows with **≥3 overlapping curated genes**, **numeric** `pathway_rank_metric` and "
            "`average_signed_effect`, AND (**pathway_rank_metric ≥ 0.5** OR **significant_gene_count ≥ 1**). "
            "Mapped-but-unscored overlap rows are in `pathway_scores_overlap_but_unscored.csv`. "
            "Primary score-based tables: `pathway_scores_evidence_filtered.csv`, `pathway_scores.csv`.*",
            "",
            "**OSD-570 PBMC:** the processed contrast is **R+45 vs R+1** — **recovery-period timing**, not flight vs pre-flight baseline.",
            "**pathway_rank_metric** is an exploratory ranking statistic — **not clinical severity** and **not DNA lesion burden**.",
            "Where overlap counts are **3–5 genes**, pathway-level signals are **limited** by coverage.",
            "",
        ]
    )

    if pl.empty:
        lines.append(
            "_No rows met the main candidate filter (overlap ≥3 and rank/significance gate) in this run._"
        )
    else:
        top = pl.sort_values(
            "pathway_rank_metric", ascending=False, na_position="last"
        ).head(12)
        for _, r in top.iterrows():
            prm = (
                float(r["pathway_rank_metric"])
                if pd.notna(r["pathway_rank_metric"])
                else float("nan")
            )
            ase = (
                float(r["average_signed_effect"])
                if pd.notna(r["average_signed_effect"])
                else float("nan")
            )
            lines.append(
                f"- **{r['dataset_name']}** | {r['comparison']} | {r['cell_type']} | "
                f"**{r['pathway']}**: exploratory ranking metric {prm:.3f} "
                f"({r.get('evidence_tier', '')}); mean signed log2 FC {ase:.4f}."
            )

    lines.extend(
        [
            "",
            "## Gene-level observations",
            "",
            "*Rows with **1–2 overlapping genes** are **gene-level observations only** — **not** pathway-level evidence.*",
            "",
        ]
    )

    if gl.empty:
        lines.append("_No 1–2 gene overlap strata in pathway_scores for this run._")
    else:
        sg = gl.sort_values(
            "pathway_rank_metric", ascending=False, na_position="last"
        ).head(15)
        for _, r in sg.iterrows():
            lines.append(
                f"- **{r['dataset_name']}** | {r['cell_type']} | {r['pathway']}: "
                f"n={int(r['overlapping_gene_count'])} — **gene-level observation only**."
            )

    lines.extend(
        [
            "",
            "## Top genes from pathway-level signals",
            "",
            "Prefer **`outputs/tables/top_genes_from_pathway_level_signals.csv`** for presentation: genes tied to "
            "**parent pathways with ≥3 overlapping genes** only (numeric `gene_signal` and `effect_size` required). "
            "**`top_damage_response_genes.csv`** includes **all** gene-level rows (including low-overlap observations).",
            "",
        ]
    )

    if (
        mapping_loaded
        and osd569_has_pathway_overlap
        and score569 is False
        and "OSD-569_whole_blood_RNA" in datasets_used
    ):
        lines.append(
            "_Mapped OSD-569 genes were excluded from top-gene summaries because numeric effect sizes were not "
            "parsed._"
        )
        lines.append("")

    if osd569_mapping_integrity and int(osd569_mapping_integrity.get("mismatched_rows", 0) or 0) > 0:
        lines.append(
            "_**OSD-569 top-gene rows were excluded** because `original_gene_id` → gene-symbol integrity "
            "validation failed (see `outputs/logs/osd569_mapping_integrity_check.txt`). Remaining rows below "
            "use the validated export table._"
        )
        lines.append("")

    pl_g = (
        pathway_level_top_genes_df
        if pathway_level_top_genes_df is not None
        else pathway_level_top_genes_table(top_genes, load_report)
    )
    if pl_g.empty:
        lines.append("_No gene rows meeting the pathway-level parent overlap filter._")
    else:
        brief = pl_g.head(12)
        for _, r in brief.iterrows():
            gs = r.get("gene_signal")
            if pd.isna(pd.to_numeric(gs, errors="coerce")):
                continue
            lines.append(
                f"- **{r['gene']}** | {r['pathway']} | {r['dataset_name']} | {r['cell_type']} | "
                f"gene_signal={float(gs):.4f}"
            )

    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- Exploratory only; **not diagnostic**.",
            "- **RNA expression does not prove DNA lesions.**",
            "- **Small crew sample size** (e.g. Inspiration4 n=4).",
            "- Immune composition and inflammation may confound pathway readouts.",
            "- **OSD-570 R+45 vs R+1** is **recovery timing**, not direct flight perturbation vs baseline.",
            "- **OSD-569** currently uses a curated-only Ensembl→HGNC map for the project's selected pathway genes.",
            "- **OSD-569** passes mapping and numeric-loader validation, but its whole-blood pathway effects are small and do not pass the main-candidate threshold.",
            "",
            "## Next steps",
            "",
            "1. Expand OSD-569 from curated-only Ensembl→HGNC mapping to a genome-wide annotation map if broader pathway libraries are added.",
            "2. Add or locate **R+1 vs pre-flight** (or similar) contrasts when available.",
            "3. Compare **flight perturbation** vs **recovery shift** explicitly.",
            "4. Expand curated gene sets using established libraries when appropriate.",
            "",
            "---",
            "",
            "_Machine-readable outputs: `outputs/logs/initial_analysis_summary.json`, `outputs/logs/data_debug_report.txt`._",
        ]
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
