#!/usr/bin/env python3
from typing import Annotated, Any, cast
from pathlib import Path

import typer
import snapatac2._snapatac2 as internal
import scanpy as sc
import anndata as ad
import numpy as np
import scipy.sparse as sp
from tqdm.auto import tqdm

app = typer.Typer()


def aggregate_sum_by_group(
    adata: ad.AnnData,
    groupby: str,
    chunk_size: int = 10000,
) -> tuple[np.ndarray, np.ndarray]:
    groups = adata.obs[groupby].to_numpy()
    if groups.size != adata.n_obs:
        raise ValueError(
            "the length of `groupby` should equal the number of observations"
        )

    unique_groups, inverse = np.unique(groups, return_inverse=True)
    n_groups = unique_groups.size
    sums = np.zeros((n_groups, adata.n_vars), dtype=np.float32)

    for chunk, start, stop in tqdm(adata.chunked_X(chunk_size), total=adata.n_obs // chunk_size + 1):
        chunk_groups = inverse[start:stop]
        present_groups = np.unique(chunk_groups)

        if sp.issparse(chunk):
            chunk = cast(Any, chunk).tocsr()
            for g in present_groups:
                rows = np.where(chunk_groups == g)[0]
                if rows.size == 0:
                    continue
                sums[g] += (
                    np.asarray(chunk[rows].sum(axis=0))
                    .ravel()
                    .astype(np.float32, copy=False)
                )
        else:
            chunk_arr = np.asarray(chunk, dtype=np.float32)
            for g in present_groups:
                rows = chunk_groups == g
                if not np.any(rows):
                    continue
                sums[g] += chunk_arr[rows].sum(axis=0, dtype=np.float32)

    return unique_groups, sums


@app.command()
def main(
    filepath: Annotated[Path, typer.Argument(help="Path to the input h5ad file")],
    output: Annotated[Path, typer.Argument(help="Path to the output h5ad file")],
    groupby: str = typer.Option(
        ..., help="Column in adata.obs to group by for feature selection."
    ),
    n_iter: int = typer.Option(
        default=1,
        help="Number of iterations to run the feature selection process",
    ),
    n_features: int = typer.Option(
        default=50000,
        help="Number of features to select based on mean accessibility (e.g. 10000).",
    ),
    blacklist: Path | None = typer.Option(
        default=None,
        help="Path to a BED file containing genomic regions to exclude from feature selection.",
    ),
    subset_selected: bool = typer.Option(
        default=False,
        help="Whether to subset the AnnData object to only include the selected features. If False, the selected features will be added to .var['selected'] and the original feature matrix will be retained.",
    ),
    chunked: bool = typer.Option(
        default=False,
        help="Whether to process the data in chunks. This can reduce memory usage for large datasets, but may increase runtime. If False, the entire dataset will be loaded into memory.",
    ),
    chunk_size: int = typer.Option(
        default=10000,
        help="Number of rows to process at a time when computing mean accessibility for feature selection. This can be adjusted based on the available memory.",
    ),
):
    # Read the input h5ad file
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    # Iteratively select features
    n_iter = max(1, n_iter)
    blacklist_mask = (
        np.array(internal.intersect_bed(adata.var_names, str(blacklist)))
        if blacklist is not None
        else None
    )

    typer.echo(
        f"Running feature selection with n_iter={n_iter}, n_features={n_features}, "
        f"blacklist={'provided' if blacklist is not None else 'not provided'}, "
        f"subset_selected={subset_selected}, chunk_size={chunk_size}"
    )

    selected_features = np.array([], dtype=np.int64)
    for _ in range(n_iter):
        typer.echo("Aggregate data by group")
        if chunked:
            _, X = aggregate_sum_by_group(adata, groupby=groupby, chunk_size=chunk_size)
        else:
            X = sc.get.aggregate(adata, by=groupby, func="sum", axis=0).layers["sum"]
        typer.echo("Calculate library sizes")
        lib_sizes = X.sum(axis=1, keepdims=True)
        lib_sizes[lib_sizes == 0] = 1.0
        typer.echo("Calculate RPM")
        rpm = X / (lib_sizes / 1_000_000.0)
        typer.echo("Calculate variance")
        var = np.var(np.log1p(rpm), axis=0)
        typer.echo("Select top features")
        selected_features = np.argsort(var)[::-1][:n_features]

        # Apply blacklist to the result
        if blacklist_mask is not None:
            typer.echo("Applying blacklist to selected features")
            selected_features = selected_features[
                np.logical_not(blacklist_mask[selected_features])
            ]

    typer.echo(f"Selected {selected_features.size} features after {n_iter} iterations")
    result = np.zeros(adata.shape[1], dtype=bool)
    result[selected_features] = True
    adata.var["selected"] = result

    # Drop unselected features if subset_selected is True
    if subset_selected:
        typer.echo("Subsetting AnnData to selected features")
        adata = adata[:, adata.var["selected"].to_numpy()]

    # Write the output h5ad file
    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()
