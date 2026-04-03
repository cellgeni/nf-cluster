#!/usr/bin/env python3
from typing import Annotated, Literal
from pathlib import Path

import typer
import snapatac2 as snap
import anndata as ad
import numpy as np
from scipy.stats import median_abs_deviation


app = typer.Typer()

genomes = {
    "hg38": snap.genome.hg38,
    "hg19": snap.genome.hg19,
    "mm10": snap.genome.mm10,
    "mm39": snap.genome.mm39,
    "GRCh37": snap.genome.GRCh37,
    "GRCh38": snap.genome.GRCh38,
    "GRCm38": snap.genome.GRCm38,
    "GRCm39": snap.genome.GRCm39,
}


@app.command()
def main(
    filepath: Annotated[Path, typer.Argument(help="Path to the input h5ad file")],
    output: Annotated[Path, typer.Argument(help="Path to the output h5ad file")],
    genome: Literal[
        "hg38", "hg19", "mm10", "mm39", "GRCh37", "GRCh38", "GRCm38", "GRCm39"
    ] = typer.Option(
        help="Genome assembly (e.g. hg38, hg19, mm10, mm39, GRCh37, GRCh38, GRCm38, GRCm39)",
    ),
):
    # Read the input h5ad file
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    # Calculate a fraction of mitochondrial fragments per cell
    if {"n_raw_reads", "n_mitochondrial_reads"}.issubset(adata.obs.columns):
        num = adata.obs["n_mitochondrial_reads"].to_numpy(dtype=float)
        den = adata.obs["n_raw_reads"].to_numpy(dtype=float)
        adata.obs["frac_mito_reads"] = np.divide(
            num,
            den,
            out=np.full_like(num, np.nan, dtype=float),  # or 0.0 if you prefer
            where=np.isfinite(den) & (den > 0),
        )

    # Calculate TSS enrichment per cell
    snap.metrics.tsse(adata, genomes[genome])

    # Calculate log1p of the number of fragments per cell
    if "n_fragments" in adata.obs.columns:
        adata.obs["log1p_n_fragments"] = np.log1p(adata.obs["n_fragments"])

    # Calculate z-score like statistic for selected metrics per cell
    columns = ["frac_mito_reads", "tss_enrichment", "log1p_n_fragments"]
    for column in columns:
        if column in adata.obs.columns:
            x = adata.obs[column].to_numpy()
            med = np.nanmedian(x)
            mad = median_abs_deviation(x, scale="normal", nan_policy="omit")
            num = x - med
            score = np.zeros_like(num, dtype=float)
            ok = np.isfinite(mad) and (mad > 1e-8)

            if ok:
                np.divide(num, mad, out=score, where=True)

            adata.obs[f"mad_score_{column}"] = score

    # Write the output h5ad file
    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()
