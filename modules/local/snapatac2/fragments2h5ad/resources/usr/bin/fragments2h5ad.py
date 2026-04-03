#!/usr/bin/env python3
import gzip
from typing import Annotated, Literal
from pathlib import Path

import typer
import snapatac2

app = typer.Typer()

genomes = {
    "hg38": snapatac2.genome.hg38,
    "hg19": snapatac2.genome.hg19,
    "mm10": snapatac2.genome.mm10,
    "mm39": snapatac2.genome.mm39,
    "GRCh37": snapatac2.genome.GRCh37,
    "GRCh38": snapatac2.genome.GRCh38,
    "GRCm38": snapatac2.genome.GRCm38,
    "GRCm39": snapatac2.genome.GRCm39,
}


@app.command()
def main(
    input: Annotated[Path, typer.Argument(help="Path to the fragment file")],
    output: Annotated[Path, typer.Argument(help="Path to the output h5ad file")],
    genome: Literal[
        "hg38", "hg19", "mm10", "mm39", "GRCh37", "GRCh38", "GRCm38", "GRCm39"
    ] = typer.Option(
        help="Genome assembly (e.g. hg38, hg19, mm10, mm39, GRCh37, GRCh38, GRCm38, GRCm39)",
    ),
    whitelist: Path = typer.Option(
        default=None,
        help="Path to a whitelist of cell barcodes. If provided, only fragments with barcodes in the whitelist will be retained. The whitelist should be a text file with one barcode per line. Can be gzipped (e.g. barcodes.tsv.gz)",
    ),
    min_num_fragments: int = typer.Option(
        default=100, help="Minimum number of fragments per cell"
    ),
    sorted_by_barcode: bool = typer.Option(
        default=False, help="Whether the fragment file is sorted by barcode"
    ),
    is_paired: bool = typer.Option(
        default=True,
        help="Indicate whether the fragment file contains paired-end reads.",
    ),
    chr_mito: Path = typer.Option(
        default=None,
        help="A list of chromosome names that are considered mitochondrial DNA. This is used to compute the fraction of mitochondrial DNA",
    ),
    sample_key: str = typer.Option(
        default=None,
        help="If specified, the sample key added to .obs.",
    ),
    n_jobs: int = typer.Option(default=-1, help="Number of parallel jobs to use"),
):
    # Get the genome object from the string
    genome = genomes.get(genome, None)
    if genome is None:
        raise typer.BadParameter("Invalid genome assembly")

    # Read mitochondrial chromosome names
    typer.echo("Reading mitochondrial chromosome names...")
    chr_mito_list = ["chrM", "M", "ChrM"]
    if chr_mito is not None:
        with open(chr_mito, "r") as f:
            chr_mito_list = [line.strip() for line in f]

    # Read the whitelist of barcodes if it is gzipped
    whitelist_barcodes = None
    if whitelist is not None and whitelist.suffix == ".gz":
        typer.echo("Reading whitelist of barcodes from gzipped file...")
        with gzip.open(whitelist, "rt") as f:
            whitelist_barcodes = set(line.strip() for line in f)

    typer.echo(
        f"Converting fragments from {input} to {output} with minimum {min_num_fragments} fragments per cell"
    )
    adata = snapatac2.pp.import_fragments(
        input,
        genome,
        whitelist=whitelist_barcodes if whitelist_barcodes is not None else whitelist,
        min_num_fragments=min_num_fragments,
        sorted_by_barcode=sorted_by_barcode,
        is_paired=is_paired,
        chrM=chr_mito_list,
        n_jobs=n_jobs,
    )

    # Add sample key to .obs if specified
    if sample_key is not None:
        adata.obs["sample"] = sample_key

    typer.echo(f"Fragments from {input} have been converted to {output}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()
