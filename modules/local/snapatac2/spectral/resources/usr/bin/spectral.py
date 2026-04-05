#!/usr/bin/env python3
from pathlib import Path
from typing import Annotated, Literal

import anndata as ad
import snapatac2 as snap
import typer


app = typer.Typer()


def _parse_sample_size(value: str | None) -> int | float | None:
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return float(value)


@app.command()
def main(
    filepath: Annotated[Path, typer.Argument(help="Path to the input h5ad file")],
    output: Annotated[Path, typer.Argument(help="Path to the output h5ad file")],
    n_comps: int = typer.Option(
        default=30,
        help="Number of dimensions to keep.",
    ),
    features: str = typer.Option(
        default="selected",
        help="Feature mask name in adata.var, or None to use all features.",
    ),
    random_state: int = typer.Option(
        default=0,
        help="Seed of the random state generator.",
    ),
    sample_size: str | None = typer.Option(
        default=None,
        help="Sample size for Nystrom approximation; integer count or fraction in (0, 1].",
    ),
    sample_method: Literal["random", "degree"] = typer.Option(
        default="random",
        help="Sampling method for Nystrom approximation.",
    ),
    chunk_size: int = typer.Option(
        default=5000,
        help="Chunk size used in the Nystrom method.",
    ),
    distance_metric: Literal["jaccard", "cosine"] = typer.Option(
        default="cosine",
        help="Distance metric for spectral embedding.",
    ),
    weighted_by_sd: bool = typer.Option(
        default=True,
        help="Whether to weight eigenvectors by square root of eigenvalues.",
    ),
    feature_weights: list[float] = typer.Option(
        default=[],
        help="Optional feature weights; repeat --feature-weights to pass multiple values.",
    ),
    num_threads: int = typer.Option(
        default=32,
        help="Number of threads to use.",
    ),
):
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    parsed_sample_size = _parse_sample_size(sample_size)
    parsed_features = None if features.lower() == "none" else features
    spectral_kwargs = {"num_threads": num_threads}

    snap.tl.spectral(
        adata,
        n_comps=n_comps,
        features=parsed_features,
        random_state=random_state,
        sample_size=parsed_sample_size,
        sample_method=sample_method,
        chunk_size=chunk_size,
        distance_metric=distance_metric,
        weighted_by_sd=weighted_by_sd,
        feature_weights=feature_weights if feature_weights else None,
        inplace=True,
        **spectral_kwargs,
    )

    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()
