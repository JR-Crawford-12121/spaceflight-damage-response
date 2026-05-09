#!/usr/bin/env python3
"""
Build follow-up priority artifacts from existing pipeline outputs (read-only).

Writes:
  outputs/logs/follow_up_priorities.md
  outputs/tables/follow_up_priorities.csv

Does not run scoring or reload raw data.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "outputs" / "tables"
LOGS = ROOT / "outputs" / "logs"

MAIN_CAND = TABLES / "main_candidate_signals.csv"
TOP_GENES = TABLES / "top_genes_from_pathway_level_signals.csv"
PATH_EV = TABLES / "pathway_scores_evidence_filtered.csv"
SUMMARY_JSON = LOGS / "initial_analysis_summary.json"

TOP_PATHWAY_N = 12
TOP_GENE_N = 12


def _safe_float(x) -> float | None:
    try:
        v = float(x)
        if math.isnan(v):
            return None
        return v
    except (TypeError, ValueError):
        return None


def _short(s: str | float | None, max_len: int = 220) -> str:
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return ""
    t = str(s).strip().replace("\n", " ")
    if len(t) > max_len:
        return t[: max_len - 3] + "..."
    return t


def _load_df(path: Path) -> pd.DataFrame:
    if not path.is_file():
        return pd.DataFrame()
    return pd.read_csv(path)


def _priority_rows_from_tables() -> pd.DataFrame:
    """Build flat rows for CSV (MD reads same source CSVs independently)."""
    mc = _load_df(MAIN_CAND)
    tg = _load_df(TOP_GENES)
    pe = _load_df(PATH_EV)

    rows: list[dict[str, str]] = []
    rank = 0

    # --- Main candidates ---
    if not mc.empty and "pathway_rank_metric" in mc.columns:
        mc2 = mc.copy()
        mc2["_prm"] = pd.to_numeric(mc2["pathway_rank_metric"], errors="coerce")
        mc2 = mc2.sort_values("_prm", ascending=False, na_position="last").head(
            TOP_PATHWAY_N
        )
        for _, r in mc2.iterrows():
            rank += 1
            prm = _safe_float(r.get("pathway_rank_metric"))
            ase = _safe_float(r.get("average_signed_effect"))
            ov = r.get("overlapping_gene_count", "")
            interp = _short(r.get("interpretation_label", ""), 160)
            caution = _short(
                r.get("evidence_tier", "") or "Exploratory; not clinical severity.",
                160,
            )
            metric = f"pathway_rank_metric={prm}; avg_signed_effect={ase}; n={ov}"
            why = (
                "High exploratory pathway rank with sufficient overlap for pathway-level context "
                f"({interp})."
            )
            rows.append(
                {
                    "priority_rank": str(rank),
                    "priority_type": "main_candidate_pathway",
                    "dataset_name": _short(r.get("dataset_name", ""), 80),
                    "comparison": _short(r.get("comparison", ""), 80),
                    "cell_type": _short(r.get("cell_type", ""), 40),
                    "pathway": _short(r.get("pathway", ""), 40),
                    "gene": "",
                    "reason_to_follow_up": _short(why, 240),
                    "supporting_metric": _short(metric, 120),
                    "caution": caution,
                }
            )

    # --- Pathway-supported genes ---
    if not tg.empty and "gene_signal" in tg.columns:
        tg2 = tg.copy()
        tg2["_gs"] = pd.to_numeric(tg2["gene_signal"], errors="coerce")
        tg2 = tg2.sort_values("_gs", ascending=False, na_position="last").head(TOP_GENE_N)
        for _, r in tg2.iterrows():
            rank += 1
            gs = _safe_float(r.get("gene_signal"))
            ef = _safe_float(r.get("effect_size"))
            metric = f"gene_signal={gs}; effect_size={ef}"
            tier = _short(r.get("evidence_tier_of_parent_pathway", ""), 120)
            why = (
                "Strong gene-level contributor under pathway overlap rules—worth closer inspection, "
                "not proof of damage."
            )
            caution = _short(
                tier + " | RNA alone cannot localize injury or prove causality.",
                240,
            )
            rows.append(
                {
                    "priority_rank": str(rank),
                    "priority_type": "pathway_supported_gene",
                    "dataset_name": _short(r.get("dataset_name", ""), 80),
                    "comparison": _short(r.get("comparison", ""), 80),
                    "cell_type": _short(r.get("cell_type", ""), 40),
                    "pathway": _short(r.get("pathway", ""), 40),
                    "gene": _short(r.get("gene", ""), 24),
                    "reason_to_follow_up": _short(why, 240),
                    "supporting_metric": _short(metric, 120),
                    "caution": caution,
                }
            )

    # --- One OSD-569 broad-context pathway row ---
    if not pe.empty:
        m569 = pe[
            pe["dataset_name"].astype(str).str.contains("OSD-569", na=False)
        ].copy()
        if not m569.empty and "pathway_rank_metric" in m569.columns:
            m569["_prm"] = pd.to_numeric(m569["pathway_rank_metric"], errors="coerce")
            r = (
                m569.sort_values("_prm", ascending=False, na_position="last")
                .iloc[0]
            )
            rank += 1
            prm = _safe_float(r.get("pathway_rank_metric"))
            ase = _safe_float(r.get("average_signed_effect"))
            ov = r.get("overlapping_gene_count", "")
            metric = f"pathway_rank_metric={prm}; avg_signed_effect={ase}; n={ov}"
            why = (
                "Whole-blood OSD-569 context: scoreable coverage with small effects—use as breadth, "
                "not the strongest candidate signal."
            )
            rows.append(
                {
                    "priority_rank": str(rank),
                    "priority_type": "osd569_broad_context",
                    "dataset_name": _short(r.get("dataset_name", ""), 80),
                    "comparison": _short(r.get("comparison", ""), 80),
                    "cell_type": _short(r.get("cell_type", ""), 40),
                    "pathway": _short(r.get("pathway", ""), 40),
                    "gene": "",
                    "reason_to_follow_up": _short(why, 240),
                    "supporting_metric": _short(metric, 120),
                    "caution": _short(
                        "Small mean effects; does not meet main-candidate threshold in typical runs.",
                        240,
                    ),
                }
            )

    # --- Validation note row (from JSON if present) ---
    integrity_note = "OSD-569 top-gene Ensembl→symbol integrity: see osd569_mapping_integrity_check.txt"
    if SUMMARY_JSON.is_file():
        try:
            data = json.loads(SUMMARY_JSON.read_text(encoding="utf-8"))
            block = data.get("osd569_mapping_integrity_status") or {}
            passed = block.get("passed")
            chk = block.get("checked_rows")
            mm = block.get("mismatched_rows")
            if passed is not None:
                integrity_note = (
                    f"osd569_mapping_integrity passed={passed}; "
                    f"checked_rows={chk}; mismatched_rows={mm}"
                )
        except (json.JSONDecodeError, OSError):
            pass

    rank += 1
    rows.append(
        {
            "priority_rank": str(rank),
            "priority_type": "validation_note",
            "dataset_name": "",
            "comparison": "",
            "cell_type": "",
            "pathway": "",
            "gene": "",
            "reason_to_follow_up": _short(
                "Row-identity checks help prevent misleading gene/ID pairings in exports.",
                240,
            ),
            "supporting_metric": _short(integrity_note, 200),
            "caution": _short(
                "Consistency checks do not validate biological mechanism.",
                200,
            ),
        }
    )

    return pd.DataFrame(rows)


def _render_markdown(mc: pd.DataFrame, tg: pd.DataFrame) -> str:
    lines: list[str] = [
        "# Follow-Up Priorities",
        "",
        "## What this report means",
        "",
        "This report does **not** identify damaged genes, DNA lesions, cancer, or clinical injury. "
        "It prioritizes RNA-expression signals in DNA damage-response pathways for **follow-up** "
        "under this repository’s exploratory scoring rules.",
        "",
        "## Highest-priority pathway signals",
        "",
        "Source: `outputs/tables/main_candidate_signals.csv` (main candidate gate applied by this pipeline).",
        "",
    ]

    if mc.empty:
        lines.append("_No main candidate rows found—run `python main.py` or check filters._")
        lines.append("")
    else:
        mc2 = mc.copy()
        mc2["_prm"] = pd.to_numeric(mc2["pathway_rank_metric"], errors="coerce")
        mc2 = mc2.sort_values("_prm", ascending=False, na_position="last").head(
            TOP_PATHWAY_N
        )
        for _, r in mc2.iterrows():
            ds = r.get("dataset_name", "")
            comp = r.get("comparison", "")
            ct = r.get("cell_type", "")
            pw = r.get("pathway", "")
            prm = r.get("pathway_rank_metric", "")
            ase = r.get("average_signed_effect", "")
            ov = r.get("overlapping_gene_count", "")
            interp = r.get("interpretation_label", "")
            tier = r.get("evidence_tier", "")
            lines.append(f"### {ds} | {ct} | {pw}")
            lines.append("")
            lines.append(f"- **Comparison:** {comp}")
            lines.append(f"- **pathway_rank_metric:** {prm}")
            lines.append(f"- **average_signed_effect:** {ase}")
            lines.append(f"- **overlapping_gene_count:** {ov}")
            lines.append(
                f"- **Why it deserves attention:** Exploratory pathway ranking with overlap ≥3 and "
                f"main-candidate filters passed; interpretation label: {interp}."
            )
            lines.append(
                f"- **Caution:** {tier or 'Exploratory metric'}; small crew sizes possible; not clinical severity."
            )
            lines.append("")

    lines.extend(
        [
            "## OSD-569 whole-blood context",
            "",
            "OSD-569 is **scoreable** in this pipeline and can show **broad** curated-gene coverage in whole blood, "
            "but observed effects are typically **small**, and **no OSD-569 pathway** passes the **main-candidate** "
            "threshold in the reference configuration. Treat OSD-569 as **broad whole-blood context**, not the "
            "strongest pathway signal on its own.",
            "",
            "## OSD-570 PBMC context",
            "",
            "OSD-570 tends to carry the **strongest candidate pathway-level signals** in this project when present, "
            "but the packaged contrast is **R+45 vs R+1**—**recovery-period timing**, not a direct "
            "flight-vs-pre-flight baseline.",
            "",
            "## Genes worth closer inspection",
            "",
            "Source: `outputs/tables/top_genes_from_pathway_level_signals.csv`. Listed genes are **worth closer "
            "inspection** under this scoring framework—they are **not** labeled as damaged genes.",
            "",
        ]
    )

    if tg.empty:
        lines.append("_No pathway-supported gene rows found._")
        lines.append("")
    else:
        tg2 = tg.copy()
        tg2["_gs"] = pd.to_numeric(tg2["gene_signal"], errors="coerce")
        tg2 = tg2.sort_values("_gs", ascending=False, na_position="last").head(TOP_GENE_N)
        for _, r in tg2.iterrows():
            gene = r.get("gene", "")
            pw = r.get("pathway", "")
            ds = r.get("dataset_name", "")
            ct = r.get("cell_type", "")
            gs = r.get("gene_signal", "")
            ef = r.get("effect_size", "")
            padj = r.get("adjusted_p_value", "")
            lines.append(f"### Gene {gene} ({pw})")
            lines.append("")
            lines.append(f"- **Dataset:** {ds} | **Cell type:** {ct}")
            lines.append(f"- **gene_signal:** {gs} | **effect_size:** {ef} | **adjusted_p_value:** {padj}")
            lines.append(
                "- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, "
                "composition); does not localize injury or predict cancer."
            )
            lines.append("")

    lines.extend(
        [
            "## Practical use",
            "",
            "Frame this project as **triage for astronaut omics follow-up**: it can help decide which **pathways**, "
            "**genes**, **immune-cell contexts**, or **assays** deserve closer attention after spaceflight.",
            "",
            "Possible follow-up assays include (examples only): **direct sequencing**, **DNA lesion assays**, "
            "**comet assay**, **gamma-H2AX staining**, **cytogenetics**, **protein-level validation**, and other "
            "orthogonal measures—not inferred from RNA alone.",
            "",
            "These signatures may help **prioritize long-term monitoring questions**, but they **do not diagnose "
            "cancer** or **predict tissue-specific disease**.",
            "",
            "## What not to overclaim",
            "",
            "- Does **not** detect DNA damage directly.",
            "- Does **not** identify damaged genes.",
            "- Does **not** predict cancer.",
            "- Does **not** localize injury to a body part.",
            "- RNA expression can reflect **inflammation**, **stress**, **immune shifts**, or **cell-composition** "
            "changes—not only DNA repair biology.",
            "",
            "---",
            "",
            "_Generated by `scripts/build_follow_up_priorities.py` from existing outputs._",
        ]
    )

    return "\n".join(lines)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    mc = _load_df(MAIN_CAND)
    tg = _load_df(TOP_GENES)

    df_out = _priority_rows_from_tables()
    if df_out.empty:
        df_out = pd.DataFrame(
            columns=[
                "priority_rank",
                "priority_type",
                "dataset_name",
                "comparison",
                "cell_type",
                "pathway",
                "gene",
                "reason_to_follow_up",
                "supporting_metric",
                "caution",
            ]
        )

    out_csv = TABLES / "follow_up_priorities.csv"
    df_out.to_csv(out_csv, index=False)

    md = _render_markdown(mc, tg)
    out_md = LOGS / "follow_up_priorities.md"
    out_md.write_text(md, encoding="utf-8")

    print(f"Wrote {out_csv.relative_to(ROOT)}")
    print(f"Wrote {out_md.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
