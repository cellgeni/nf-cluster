#!/usr/bin/env python3
from typing import Annotated, Literal
from pathlib import Path

import typer
import snapatac2 as snap
import anndata as ad


app = typer.Typer()


@app.command()
def main(
    filepath: Annotated[Path, typer.Argument(help="Path to the input h5ad file")],
    output: Annotated[Path, typer.Argument(help="Path to the output h5ad file")],
    bin_size: int = typer.Option(
        default=500,
        help="Bin size for tiling the genome (e.g. 500 for 500bp bins).",
    ),
    exclude_chroms: list[str] = typer.Option(
        default=["chrM", "chrY", "M", "Y"],
        help="List of chromosome names to exclude when adding tile information (e.g. ['chrM', 'chrY', 'M', 'Y']).",
    ),
    min_frag_size: int = typer.Option(
        default=None,
        help="Minimum fragment size for filtering (e.g. 50 for 50bp fragments).",
    ),
    max_frag_size: int = typer.Option(
        default=None,
        help="Maximum fragment size for filtering (e.g. 1000 for 1000bp fragments).",
    ),
    counting_strategy: Literal[
        "fragment", "insertion", "paired-insertion"
    ] = typer.Option(
        default="paired-insertion",
        help="Counting strategy for adding tile information (e.g. 'fragment', 'insertion', 'paired-insertion').",
    ),
    value_type: Literal["target", "total", "fraction"] = typer.Option(
        default="target",
        help="The type of value to use from .obsm['_values']",
    ),
    summary_type: Literal["sum", "mean"] = typer.Option(
        default="sum",
        help="Summary type to use when counting_strategy is 'paired-insertion' (e.g. 'sum' or 'mean').",
    ),
    
):
    # Read the input h5ad file
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    # Add tile information to the AnnData object
    snap.pp.add_tile_matrix(
        adata,
        bin_size=bin_size,
        exclude_chroms=exclude_chroms,
        min_frag_size=min_frag_size,
        max_frag_size=max_frag_size,
        counting_strategy=counting_strategy,
        summary_type=summary_type,
        value_type=value_type,
        inplace=True,
    )

    # Write the output h5ad file
    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()
