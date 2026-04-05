#!/usr/bin/env python3
from pathlib import Path
from typing import Annotated, Literal
import importlib

import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
import cupy as cp
import anndata as ad
import rapids_singlecell as rsc
import typer

rmm.initialize()
cp.cuda.set_allocator(rmm_cupy_allocator)

app = typer.Typer()


@app.command()
def main(
    filepath: Annotated[Path, typer.Argument(help="Path to the input h5ad file")],
    output: Annotated[Path, typer.Argument(help="Path to the output h5ad file")],
    n_neighbors: int = typer.Option(
        default=15,
        help="Size of local neighborhood.",
    ),
    n_pcs: int | None = typer.Option(
        default=None,
        help="Number of principal components to use. If None, defaults are used by rapids-singlecell.",
    ),
    use_rep: str | None = typer.Option(
        default=None,
        help="Representation key to use (e.g. X_spectral).",
    ),
    random_state: int = typer.Option(
        default=0,
        help="Random state for neighbor graph computation.",
    ),
    algorithm: Literal[
        "brute",
        "cagra",
        "ivfflat",
        "ivfpq",
        "nn_descent",
        "all_neighbors",
        "mg_ivfflat",
        "mg_ivfpq",
    ] = typer.Option(
        default="brute",
        help="Neighbor search algorithm.",
    ),
    metric: str = typer.Option(
        default="euclidean",
        help="Distance metric.",
    ),
    method: Literal["umap", "gauss", "jaccard"] = typer.Option(
        default="umap",
        help="Method for computing connectivities.",
    ),
    key_added: str | None = typer.Option(
        default=None,
        help="Optional key under .uns where neighbors metadata is stored.",
    ),
):
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    # Transfer data to GPU
    rsc.get.anndata_to_GPU(adata)

    # Compute neighbors using rapids-singlecell
    rsc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep,
        random_state=random_state,
        algorithm=algorithm,
        metric=metric,
        method=method,
        key_added=key_added,
        copy=False,
    )

    # Transfer results back to CPU
    rsc.get.anndata_to_CPU(adata)

    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()
