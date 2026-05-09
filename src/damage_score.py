"""
Standardize differential expression rows and compute pathway-level scores.

This quantifies exploratory transcriptional signatures associated with DNA
damage-response biology — not direct DNA lesion detection.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from gene_sets import PATHWAY_GENE_SETS
from utils import safe_float, significance_weight


def normalize_gene_symbol(g: str) -> str:
    return str(g).strip().upper()


def evidence_tier_label(overlapping_gene_count: int) -> str:
    """Human-readable evidence strength from overlap count only."""
    n = int(overlapping_gene_count)
    if n == 0:
        return "no overlapping genes"
    if n in (1, 2):
        return "gene-level signal only"
    if 3 <= n <= 5:
        return "limited pathway-level signal"
    return "stronger pathway-level coverage"


def compute_gene_signals(df: pd.DataFrame) -> pd.Series:
    """gene_signal = abs(effect_size) * significance_weight(adjusted_p)."""
    signals = []
    for _, row in df.iterrows():
        eff = safe_float(row.get("effect_size"))
        if eff is None:
            signals.append(np.nan)
            continue
        padj = safe_float(row.get("adjusted_p_value"))
        w = significance_weight(padj if padj is not None else None)
        signals.append(abs(eff) * w)
    return pd.Series(signals, index=df.index)


def score_one_pathway(
    standardized_df: pd.DataFrame,
    pathway_name: str,
    genes_in_pathway: set[str],
) -> dict[str, Any]:
    """Aggregate metrics for one pathway gene set."""
    sub = standardized_df[
        standardized_df["gene"].map(normalize_gene_symbol).isin(genes_in_pathway)
    ].copy()

    if sub.empty:
        return {
            "pathway": pathway_name,
            "overlapping_gene_count": 0,
            "pathway_rank_metric": np.nan,
            "average_signed_effect": np.nan,
            "average_absolute_effect": np.nan,
            "upregulated_gene_count": 0,
            "downregulated_gene_count": 0,
            "significant_gene_count": 0,
            "direction": "unknown",
            "top_genes_detail": pd.DataFrame(),
        }

    sub["gene_signal"] = compute_gene_signals(sub)
    sub["direction"] = np.where(sub["effect_size"] >= 0, "up", "down")

    signed = pd.to_numeric(sub["effect_size"], errors="coerce")
    avg_signed = float(np.nanmean(signed))
    avg_abs = float(np.nanmean(np.abs(signed)))

    up_n = int((signed > 0).sum())
    down_n = int((signed < 0).sum())

    padj = pd.to_numeric(sub["adjusted_p_value"], errors="coerce")
    sig_n = int((padj < 0.05).sum()) if padj.notna().any() else 0

    pathway_rank_metric = float(np.nanmean(sub["gene_signal"]))

    top_detail = sub.sort_values("gene_signal", ascending=False).copy()
    top_detail.insert(0, "pathway", pathway_name)
    tier_parent = evidence_tier_label(len(sub))
    top_detail["evidence_tier_of_parent_pathway"] = tier_parent
    top_detail["parent_pathway_overlap_count"] = len(sub)
    # Do NOT reassign original_gene_id from unsorted `sub` after sort_values — that scrambles
    # row identity vs gene / effect_size (silent bug). Rows stay aligned from sorted top_detail.
    if "original_gene_id" not in top_detail.columns:
        top_detail["original_gene_id"] = top_detail["gene"].astype(str)

    return {
        "pathway": pathway_name,
        "overlapping_gene_count": len(sub),
        "pathway_rank_metric": pathway_rank_metric,
        "average_signed_effect": avg_signed,
        "average_absolute_effect": avg_abs,
        "upregulated_gene_count": up_n,
        "downregulated_gene_count": down_n,
        "significant_gene_count": sig_n,
        "direction": "positive" if avg_signed > 0 else "negative",
        "top_genes_detail": top_detail,
    }


def score_all_pathways(standardized_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      pathway_scores_df — one row per pathway per (dataset_name, comparison, cell_type)
      top_genes_df — ranked gene-level contributors
    """
    pathways_normalized = {
        name: {normalize_gene_symbol(g) for g in genes}
        for name, genes in PATHWAY_GENE_SETS.items()
    }

    # Stratify so PBMC cell types stay separate when present
    group_cols = ["dataset_name", "comparison", "cell_type"]
    rows_scores: list[dict[str, Any]] = []
    top_chunks: list[pd.DataFrame] = []

    for (ds, comp, ct), grp in standardized_df.groupby(group_cols, dropna=False):
        for pname, gset in pathways_normalized.items():
            res = score_one_pathway(grp, pname, gset)
            detail = res.pop("top_genes_detail")
            res["dataset_name"] = ds
            res["comparison"] = comp
            res["cell_type"] = ct if pd.notna(ct) else "unknown"

            # Flatten top genes for export
            if not detail.empty:
                top_chunks.append(detail)

            rows_scores.append(res)

    pathway_scores_df = pd.DataFrame(rows_scores)

    keep_cols = [
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

    if top_chunks:
        top_genes_df = pd.concat(top_chunks, ignore_index=True)
        for c in keep_cols:
            if c not in top_genes_df.columns:
                top_genes_df[c] = np.nan
        top_genes_df = top_genes_df[keep_cols]
        top_genes_df = top_genes_df.sort_values(
            ["pathway", "gene_signal"], ascending=[True, False]
        )
    else:
        top_genes_df = pd.DataFrame(columns=keep_cols)

    return pathway_scores_df, top_genes_df


# --- Interpretation thresholds (edit freely) ---
MIN_OVERLAP_GENES = 3
STRONG_POS_SCORE = 2.0
MODERATE_POS_SCORE = 1.0


def interpretation_label(row: pd.Series) -> str:
    """Conservative text labels for pathway_rank_metric rows."""
    n = int(row.get("overlapping_gene_count", 0) or 0)
    if n < MIN_OVERLAP_GENES:
        return "insufficient overlapping genes"

    score = safe_float(row.get("pathway_rank_metric"))
    avg_eff = safe_float(row.get("average_signed_effect"))
    if score is None or avg_eff is None:
        return "weak or unclear signature"

    if score >= STRONG_POS_SCORE and avg_eff > 0:
        return "strong positive DNA damage-response signature"
    if score >= MODERATE_POS_SCORE and avg_eff > 0:
        return "moderate positive DNA damage-response signature"
    if score >= MODERATE_POS_SCORE and avg_eff < 0:
        return "pathway genes changed but mostly downregulated"
    return "weak or unclear signature"


def attach_interpretation_labels(pathway_scores_df: pd.DataFrame) -> pd.DataFrame:
    out = pathway_scores_df.copy()
    out["interpretation_label"] = out.apply(interpretation_label, axis=1)
    out["evidence_tier"] = out["overlapping_gene_count"].apply(evidence_tier_label)
    return out
