"""
Matplotlib / seaborn figures for pathway overview (publication-style simplicity).

Primary pathway plots use overlap ≥3 only and exclude datasets with no curated-symbol overlap.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from utils import project_root

MIN_OVERLAP_FOR_PRIMARY_PLOTS = 3


def _fig_path(name: str, root: Path | None = None) -> Path:
    root = root or project_root()
    return root / "outputs" / "figures" / name


def _plot_ready_pathway_rows(
    pathway_scores: pd.DataFrame,
    excluded_plot_datasets: set[str],
) -> pd.DataFrame:
    """Rows eligible for main pathway figures (not severity charts)."""
    df = pathway_scores.copy()
    df = df[df["overlapping_gene_count"] >= MIN_OVERLAP_FOR_PRIMARY_PLOTS]
    if excluded_plot_datasets:
        df = df[~df["dataset_name"].astype(str).isin(excluded_plot_datasets)]
    return df


def plot_pathway_scores_bar(
    pathway_scores: pd.DataFrame,
    excluded_plot_datasets: set[str],
    root: Path | None = None,
) -> None:
    """Bar chart of exploratory ranking metric by pathway (overlap ≥3 only)."""
    root = root or project_root()
    df = _plot_ready_pathway_rows(pathway_scores, excluded_plot_datasets)
    if df.empty:
        return

    df = df.copy()
    groups = list(df.groupby(["dataset_name", "comparison"], sort=False))
    n_panels = len(groups)
    fig, axes = plt.subplots(
        n_panels, 1, figsize=(11, 3.6 * max(1, n_panels)), sharex=False
    )
    if n_panels == 1:
        axes = np.array([axes])
    for ax, ((ds, comp), sub) in zip(axes, groups):
        order = sorted(sub["pathway"].unique())
        sns.barplot(
            data=sub,
            x="pathway",
            y="pathway_rank_metric",
            hue="cell_type",
            order=order,
            ax=ax,
        )
        ax.tick_params(axis="x", rotation=25)
        for lab in ax.get_xticklabels():
            lab.set_ha("right")
        ax.set_ylabel("Exploratory pathway ranking metric")
        ax.set_xlabel("Pathway")
        ax.set_title(
            f"{ds}\nContrast: {comp}\n"
            f"(pathways with ≥{MIN_OVERLAP_FOR_PRIMARY_PLOTS} overlapping genes only)"
        )
        if sub["cell_type"].nunique() <= 1:
            leg = ax.get_legend()
            if leg is not None:
                leg.remove()
        else:
            ax.legend(title="cell_type", fontsize=7, title_fontsize=8, loc="upper right")
    fig.suptitle(
        "Exploratory DNA damage-response pathway ranking\n"
        "(not clinical severity)",
        fontsize=12,
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(
        _fig_path("pathway_scores_bar_chart.png", root), dpi=150, bbox_inches="tight"
    )
    plt.close()


def plot_top_genes_heatmap(
    top_genes: pd.DataFrame,
    excluded_plot_datasets: set[str],
    root: Path | None = None,
) -> None:
    """Heatmap of effect sizes for top genes."""
    root = root or project_root()
    if top_genes.empty or "effect_size" not in top_genes.columns:
        return

    tg = top_genes.copy()
    if excluded_plot_datasets and "dataset_name" in tg.columns:
        tg = tg[~tg["dataset_name"].astype(str).isin(excluded_plot_datasets)]
    if tg.empty:
        return

    tg["_abs"] = tg["effect_size"].abs()
    top_symbols = (
        tg.groupby("gene")["_abs"].max().sort_values(ascending=False).head(20).index
    )
    sub = tg[tg["gene"].isin(top_symbols)]
    if sub.empty:
        return

    pivot = sub.pivot_table(
        index="gene",
        columns="pathway",
        values="effect_size",
        aggfunc="mean",
    )
    pivot = pivot.reindex(index=list(top_symbols))

    plt.figure(figsize=(10, max(4, 0.35 * len(pivot))))
    sns.heatmap(pivot, cmap="RdBu_r", center=0, linewidths=0.5)
    plt.title("Effect sizes for top DNA damage-response genes (mean by pathway)")
    plt.tight_layout()
    plt.savefig(_fig_path("top_genes_heatmap.png", root), dpi=150)
    plt.close()


def plot_pathway_direction(
    pathway_scores: pd.DataFrame,
    excluded_plot_datasets: set[str],
    root: Path | None = None,
) -> None:
    """Average signed effect by pathway (overlap ≥3 only)."""
    root = root or project_root()
    df = _plot_ready_pathway_rows(pathway_scores, excluded_plot_datasets)
    if df.empty:
        return

    plt.figure(figsize=(9, 4))
    order = sorted(df["pathway"].unique())
    sns.barplot(
        data=df,
        x="pathway",
        y="average_signed_effect",
        hue="dataset_name",
        order=order,
    )
    plt.axhline(0, color="gray", linewidth=0.8)
    plt.xticks(rotation=25, ha="right")
    plt.ylabel("Mean signed effect (log2 fold change)")
    plt.title(
        "Direction summary (overlap ≥"
        f"{MIN_OVERLAP_FOR_PRIMARY_PLOTS} genes only)"
    )
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(_fig_path("pathway_direction_summary.png", root), dpi=150)
    plt.close()


def plot_cell_type_pathway_scores(
    pathway_scores: pd.DataFrame,
    excluded_plot_datasets: set[str],
    root: Path | None = None,
) -> None:
    """Optional: ranking metric across immune cell types (overlap ≥3)."""
    root = root or project_root()
    df = _plot_ready_pathway_rows(pathway_scores, excluded_plot_datasets)
    df = df[df["cell_type"].astype(str) != "unknown"]
    if df.empty or df["cell_type"].nunique() < 2:
        return

    plt.figure(figsize=(10, 5))
    sns.barplot(
        data=df,
        x="pathway",
        y="pathway_rank_metric",
        hue="cell_type",
    )
    plt.xticks(rotation=25, ha="right")
    plt.ylabel("Exploratory pathway ranking metric")
    plt.title(
        "Cell-type ranking metrics "
        f"(overlap ≥{MIN_OVERLAP_FOR_PRIMARY_PLOTS} genes)"
    )
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(_fig_path("cell_type_pathway_scores.png", root), dpi=150)
    plt.close()


def plot_evidence_filtered_pathway_heatmap(
    pathway_scores: pd.DataFrame,
    excluded_plot_datasets: set[str],
    root: Path | None = None,
    footer_warning: str | None = None,
) -> None:
    """
    Heatmap: rows = cell_type, columns = pathway, color = average_signed_effect,
    annotated with n=<overlapping_gene_count>. One panel per (dataset_name, comparison).
    """
    root = root or project_root()
    df = pathway_scores.copy()
    df = df[df["overlapping_gene_count"] >= MIN_OVERLAP_FOR_PRIMARY_PLOTS]
    if excluded_plot_datasets:
        df = df[~df["dataset_name"].astype(str).isin(excluded_plot_datasets)]
    if df.empty:
        return

    groups = list(df.groupby(["dataset_name", "comparison"], sort=False))
    n_panels = len(groups)
    # Shared symmetric scale across panels so small effects (e.g. whole blood) are not
    # over-saturated relative to larger contrasts in another panel.
    max_abs = 0.0
    for _, sub in groups:
        v = pd.to_numeric(sub["average_signed_effect"], errors="coerce").to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        if v.size:
            max_abs = max(max_abs, float(np.nanmax(np.abs(v))))
    if max_abs == 0.0:
        max_abs = 1e-9
    vmin, vmax = -max_abs, max_abs

    fig, axes = plt.subplots(n_panels, 1, figsize=(11, 4 * max(1, n_panels)))
    if n_panels == 1:
        axes = np.array([axes])
    for panel_i, (ax, ((ds, comp), sub)) in enumerate(zip(axes, groups)):
        cts = sorted(sub["cell_type"].astype(str).unique())
        pws = sorted(sub["pathway"].unique())
        pv = sub.pivot_table(
            index="cell_type",
            columns="pathway",
            values="average_signed_effect",
            aggfunc="first",
        ).reindex(index=cts, columns=pws)
        pn = sub.pivot_table(
            index="cell_type",
            columns="pathway",
            values="overlapping_gene_count",
            aggfunc="first",
        ).reindex(index=cts, columns=pws)

        mask = pv.isna()
        annot = pn.map(lambda x: f"n={int(x)}" if pd.notna(x) else "")
        sns.heatmap(
            pv.astype(float),
            mask=mask,
            annot=annot,
            fmt="",
            cmap="RdBu_r",
            center=0,
            vmin=vmin,
            vmax=vmax,
            linewidths=0.5,
            ax=ax,
            cbar=panel_i == n_panels - 1,
            cbar_kws={"label": "Mean signed log2 fold change"},
        )
        # Panel label only when multiple datasets/contrasts need disambiguation
        if n_panels > 1:
            ax.set_title(f"{ds}\n{comp}", fontsize=9)
        ax.set_xlabel("pathway")
        ax.set_ylabel("cell_type")

    subtitle = ""
    if n_panels == 1:
        ds0, comp0 = groups[0][0]
        ds_s = str(ds0)
        comp_s = str(comp0)
        if "OSD-570" in ds_s and "R+45" in comp_s and "R+1" in comp_s:
            subtitle = (
                "OSD-570 PBMC snRNA | R+45 vs R+1 recovery contrast | "
                f"overlap ≥{MIN_OVERLAP_FOR_PRIMARY_PLOTS} curated genes"
            )
        else:
            subtitle = (
                f"{ds_s} | {comp_s} | overlap ≥{MIN_OVERLAP_FOR_PRIMARY_PLOTS} curated genes"
            )
    else:
        subtitle = (
            f"One panel per dataset and contrast | overlap ≥{MIN_OVERLAP_FOR_PRIMARY_PLOTS} curated genes"
        )

    fig.suptitle(
        "Evidence-filtered pathway heatmap\n" + subtitle,
        fontsize=11,
        y=1.02,
    )
    plt.tight_layout(rect=(0, 0.08, 1, 1))
    footer_txt = (
        "Shared color scale across panels. Color = mean signed log2 fold change. "
        "Blank cells = insufficient overlap. Exploratory only, not clinical severity."
    )
    if footer_warning:
        footer_txt = footer_txt + " " + footer_warning.strip()
    fig.text(
        0.5,
        0.02,
        footer_txt,
        ha="center",
        va="bottom",
        fontsize=9,
        transform=fig.transFigure,
    )
    plt.savefig(
        _fig_path("evidence_filtered_pathway_heatmap.png", root),
        dpi=150,
        bbox_inches="tight",
        pad_inches=0.35,
    )
    plt.close()


def run_all_plots(
    pathway_scores: pd.DataFrame,
    top_genes: pd.DataFrame,
    excluded_plot_datasets: set[str],
    root: Path | None = None,
    heatmap_extra_exclusions: set[str] | None = None,
    heatmap_footer_warning: str | None = None,
) -> None:
    """Generate figures; primary pathway plots omit overlap<3 and identifier-failed datasets."""
    root = root or project_root()
    (root / "outputs" / "figures").mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    plot_pathway_scores_bar(pathway_scores, excluded_plot_datasets, root)
    plot_top_genes_heatmap(top_genes, excluded_plot_datasets, root)
    plot_pathway_direction(pathway_scores, excluded_plot_datasets, root)
    plot_cell_type_pathway_scores(pathway_scores, excluded_plot_datasets, root)
    hm_exc = set(excluded_plot_datasets or [])
    if heatmap_extra_exclusions:
        hm_exc |= heatmap_extra_exclusions
    plot_evidence_filtered_pathway_heatmap(
        pathway_scores,
        hm_exc,
        root,
        footer_warning=heatmap_footer_warning,
    )
