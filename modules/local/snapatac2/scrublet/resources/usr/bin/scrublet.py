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
    features: str = typer.Option(
        default="selected",
        help="Feature set stored in adata.var to use for scrublet.",
    ),
    n_comps: int = typer.Option(
        default=15,
        help="Number of principal components used by scrublet.",
    ),
    sim_doublet_ratio: float = typer.Option(
        default=2.0,
        help="Ratio of simulated doublets to observed cells.",
    ),
    expected_doublet_rate: float = typer.Option(
        default=0.1,
        help="Expected doublet rate.",
    ),
    n_neighbors: int | None = typer.Option(
        default=None,
        help="Number of neighbors for scrublet. If unset, library default is used.",
    ),
    use_approx_neighbors: bool = typer.Option(
        default=False,
        help="Whether to use approximate nearest-neighbor search.",
    ),
    random_state: int = typer.Option(
        default=4,
        help="Random seed for scrublet.",
    ),
    probability_threshold: float = typer.Option(
        default=0.5,
        help="Probability threshold for filter_doublets.",
    ),
    score_threshold: float | None = typer.Option(
        default=None,
        help="Optional score threshold for filter_doublets.",
    ),
    filter_inplace: bool = typer.Option(
        default=False,
        help="Whether to filter doublets from the AnnData object. If False, doublet annotations will be added to adata.obs but no cells will be removed.",
    ),
):
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    # Run scrublet to identify doublets
    snap.pp.scrublet(
        adata,
        features=features,
        n_comps=n_comps,
        sim_doublet_ratio=sim_doublet_ratio,
        expected_doublet_rate=expected_doublet_rate,
        n_neighbors=n_neighbors,
        use_approx_neighbors=use_approx_neighbors,
        random_state=random_state,
    )

    # Filter doublets based on the scrublet results
    if filter_inplace:
        snap.pp.filter_doublets(
            adata,
            probability_threshold=probability_threshold,
            score_threshold=score_threshold,
            inplace=filter_inplace,
        )
    else:
        adata.obs["doublet_filter"] = snap.pp.filter_doublets(
            adata,
            probability_threshold=probability_threshold,
            score_threshold=score_threshold,
            inplace=filter_inplace,
        )

    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()
