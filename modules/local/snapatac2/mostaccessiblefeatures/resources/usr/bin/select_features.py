#!/usr/bin/env python3
from typing import Annotated
from pathlib import Path

import typer
import numpy as np
import anndata as ad
import snapatac2._snapatac2 as internal


app = typer.Typer(add_completion=False)


def _find_most_accessible_features(
    feature_count: np.ndarray,
    filter_lower_quantile: float,
    filter_upper_quantile: float,
    total_features: int,
) -> np.ndarray:
    idx = np.argsort(feature_count)

    for i in range(idx.size):
        if feature_count[idx[i]] > 0:
            break
    else:
        return np.array([], dtype=int)

    idx = idx[i:]
    n = idx.size
    n_lower = int(filter_lower_quantile * n)
    n_upper = int(filter_upper_quantile * n)
    idx = idx[n_lower : n - n_upper]

    return idx[::-1][:total_features]


def select_most_accessible_features(
    adata: ad.AnnData,
    n_features: int,
    filter_lower_quantile: float,
    filter_upper_quantile: float,
    blacklist: Path | None = None,
) -> np.ndarray:
    count = np.zeros(adata.shape[1], dtype=float)

    for batch, _, _ in adata.chunked_X(2000):
        count += np.ravel(batch.sum(axis=0))

    adata.var["count"] = count

    selected_features = _find_most_accessible_features(
        count,
        filter_lower_quantile,
        filter_upper_quantile,
        n_features,
    )

    if blacklist is not None:
        blacklist_mask = np.array(
            internal.intersect_bed(adata.var_names, str(blacklist))
        )
        selected_features = selected_features[
            np.logical_not(blacklist_mask[selected_features])
        ]

    selected_mask = np.zeros(adata.shape[1], dtype=bool)
    selected_mask[selected_features] = True
    adata.var["selected"] = selected_mask

    return selected_mask


@app.command()
def main(
    filepath: Annotated[Path, typer.Argument(help="Path to the input h5ad file")],
    output: Annotated[Path, typer.Argument(help="Path to the output h5ad file")],
    n_features: int = typer.Option(
        default=50000,
        help="Number of most accessible features to keep after filtering.",
    ),
    filter_lower_quantile: float = typer.Option(
        default=0.005,
        help="Lower quantile of nonzero feature counts to discard.",
    ),
    filter_upper_quantile: float = typer.Option(
        default=0.005,
        help="Upper quantile of nonzero feature counts to discard.",
    ),
    blacklist: Path | None = typer.Option(
        default=None,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Path to a BED file containing genomic regions to exclude.",
    ),
    subset_selected: bool = typer.Option(
        default=False,
        help="Whether to subset the AnnData object to selected features only.",
    ),
):
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    selected_mask = select_most_accessible_features(
        adata=adata,
        n_features=n_features,
        filter_lower_quantile=filter_lower_quantile,
        filter_upper_quantile=filter_upper_quantile,
        blacklist=blacklist,
    )

    typer.echo(f"Selected {int(selected_mask.sum())} features.")

    if subset_selected:
        adata = adata[:, selected_mask].copy()
        typer.echo("Subsetted AnnData object to selected features only.")

    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()