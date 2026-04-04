#!/usr/bin/env python3
import os
from typing import Annotated, Literal
from pathlib import Path

import typer
import snapatac2 as snap
import anndata as ad
import numpy as np
from scipy.stats import median_abs_deviation
import matplotlib.pyplot as plt
import seaborn as sns

app = typer.Typer()


def mad_score(x: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    """
    Calculate the median absolute deviation (MAD) score for an array of values.
    Args:
        x (np.ndarray): Array of values for which to calculate the MAD score.
        mask (np.ndarray, optional): Boolean array indicating which values to include in the MAD calculation. If None, all values are included. Defaults to None.

    Returns:
        np.ndarray: Array of MAD scores.
    """
    if mask is None:
        mask = np.isfinite(x)
    med = np.nanmedian(x[mask])
    mad = median_abs_deviation(x[mask], scale="normal", nan_policy="omit")
    num = x - med
    score = np.zeros_like(num, dtype=float)
    ok = np.isfinite(mad) and (mad > 1e-8)

    if ok:
        np.divide(num, mad, out=score, where=True)
    return score


@app.command()
def main(
    filepath: Annotated[Path, typer.Argument(help="Path to the input h5ad file")],
    output_prefix: Annotated[Path, typer.Argument(help="Prefix for the output files")],
    mad_thresh: float = typer.Option(
        default=3.5,
        help="Median absolute deviation (MAD) threshold for filtering cells based on selected metrics. Cells with a MAD value greater than mad_thresh for any of the selected metrics will be filtered out. Set to 0 to disable MAD-based filtering.",
    ),
    amulet_qthresh: float = typer.Option(
        default=0.01,
        help="Q-value threshold for AMULET doublet detection. Cells with amulet_qvalue < amulet_qthresh will be filtered out. Set to 0 to disable AMULET doublet filtering.",
    ),
    min_fragments: int = typer.Option(
        default=100,
        help="Minimum number of fragments per cell. Cells with n_fragments_cr < min_fragments will be filtered out. Set to 0 to disable filtering based on the number of fragments.",
    ),
    max_mito_frac: float = typer.Option(
        default=0.2,
        help="Maximum fraction of mitochondrial reads per cell. Cells with frac_mito_reads > max_mito_frac will be filtered out. Set to 1 to disable filtering based on the fraction of mitochondrial reads.",
    ),
    mad_key: list[str] = typer.Option(
        default=["log1p_n_fragments_cr", "frac_mito_reads", "tsse"],
        help="List of metric keys to use for MAD-based filtering. The specified keys should be present in adata.obs.",
    ),
    filter_cells: bool = typer.Option(
        default=False,
        help="Whether to filter cells based on the specified criteria. If False, the filtering results will be added to the AnnData object, but no cells will be removed.",
    ),
):
    # Read the input h5ad file
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    # Filter AMULET doublets
    failed_cells = dict()
    amulet_mask = np.full(adata.n_obs, True, dtype=bool)
    if "amulet_qvalue" in adata.obs.columns:
        amulet_doublets = adata.obs["amulet_qvalue"] < amulet_qthresh
        failed_cells["amulet"] = adata.obs[amulet_doublets].index.tolist()
        amulet_mask = (~amulet_doublets).to_numpy()

    # Calculate MAD scores and filter cells based on the specified metrics
    for key in mad_key:
        if key in adata.obs.columns:
            adata.obs[f"mad_{key}"] = mad_score(adata.obs[key].to_numpy(), amulet_mask)
            failed_cells[key] = adata.obs[
                adata.obs[f"mad_{key}"] > mad_thresh
            ].index.tolist()

    # Visualize the filtering results
    vis_cat = list(filter(lambda x: x != "amulet", failed_cells.keys()))
    tsee_frag_plot = int({"tsse", "log1p_n_fragments_cr"}.issubset(adata.obs.columns))
    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(vis_cat) + tsee_frag_plot,
        figsize=(5 * (len(vis_cat) + tsee_frag_plot), 5),
    )
    for ax, cat in zip(axes.flatten(), vis_cat):
        visdf = adata.obs[[cat]].copy()
        visdf.loc[:, "qcategory"] = "passed"
        visdf.loc[failed_cells[cat], "qcategory"] = "failed"
        sns.histplot(
            visdf,
            x=cat,
            hue="qcategory",
            ax=ax,
            palette={"failed": "red", "passed": "tab:blue"},
        )
        ax.set_title(output_prefix)
        ax.legend_.remove()

    # Add filtering results to AnnData object
    adata.obs["quality_control"] = "pass"
    for cat, cells in failed_cells.items():
        adata.obs.loc[cells, "quality_control"] = "fail"
        typer.echo(f"Filtering category: {cat}, failed cells: {len(cells)}")
    adata.uns["qc_summary"] = {cat: len(cells) for cat, cells in failed_cells.items()}
    adata.uns["qc_summary"]["total_filtered"] = (
        adata.obs["quality_control"] == "fail"
    ).sum()
    typer.echo(
        f"Total filtered cells: {adata.uns['qc_summary']['total_filtered']}/{adata.n_obs}"
    )

    # Add tsee vs n_fragments histogram
    if tsee_frag_plot:
        sns.histplot(
            adata.obs,
            x="log1p_n_fragments_cr",
            y="tsse",
            hue="quality_control",
            palette={"fail": "tab:red", "pass": "tab:blue"},
            ax=axes.flatten()[-1],
        )
        sns.kdeplot(
            adata.obs,
            x="log1p_n_fragments_cr",
            y="tsse",
            ax=axes.flatten()[-1],
        )
        axes.flatten()[-1].set_title(output_prefix)
        handles, labels = axes.flatten()[-1].get_legend_handles_labels()
        if handles:
            axes.flatten()[-1].legend(
                handles,
                labels,
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                borderaxespad=0,
            )

    # Save the filtering results visualization    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", bbox_inches="tight")

    # Save the filtered AnnData object
    if filter_cells:
        adata = adata[adata.obs["quality_control"] == "pass"].copy()
    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(f"{output_prefix}_qc.h5ad")


if __name__ == "__main__":
    app()
