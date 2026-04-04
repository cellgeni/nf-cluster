#!/usr/bin/env python3
from typing import Annotated
from pathlib import Path

import typer
import snapatac2 as snap
import anndata as ad


app = typer.Typer()


@app.command()
def main(
    filepath: Annotated[Path, typer.Argument(help="Path to the input h5ad file")],
    output: Annotated[Path, typer.Argument(help="Path to the output h5ad file")],
    n_features: int = typer.Option(
        default=50000,
        help="Number of features to select based on mean accessibility (e.g. 10000).",
    ),
    filter_lower_quantile: float = typer.Option(
        default=0.005,
        help="Lower quantile for filtering features based on mean accessibility.",
    ),
    filter_upper_quantile: float = typer.Option(
        default=0.005,
        help="Upper quantile of the feature count distribution to filter out.",
    ),
    blacklist: Path | None = typer.Option(
        default=None,
        help="Path to a BED file containing genomic regions to exclude from feature selection.",
    ),
    subset_selected: bool = typer.Option(
        default=False,
        help="Whether to subset the AnnData object to only include the selected features. If False, the selected features will be added to .var['selected'] and the original feature matrix will be retained.",
    ),
):
    # Read the input h5ad file
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    # Select features from the existing tile matrix
    snap.pp.select_features(
        adata,
        n_features=n_features,
        filter_lower_quantile=filter_lower_quantile,
        filter_upper_quantile=filter_upper_quantile,
        blacklist=blacklist,
        max_iter=1,
        inplace=True,
    )

    # Drop unselected features if subset_selected is True
    if subset_selected:
        adata = adata[:, adata.var["selected"]].copy()

    # Write the output h5ad file
    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()
