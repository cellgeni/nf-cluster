#!/usr/bin/env python3
from pathlib import Path
from typing import Annotated, Literal

import anndata as ad
import typer


app = typer.Typer()


@app.command()
def main(
    in_files: Annotated[list[Path], typer.Argument(help="Input AnnData .h5ad files to concatenate")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output concatenated .h5ad file")],
    max_loaded_elems: int = typer.Option(
        default=100000000,
        help="Maximum elements loaded in memory while concatenating sparse arrays.",
    ),
    axis: Literal["obs", "var"] = typer.Option(
        default="obs",
        help="Axis to concatenate along.",
    ),
    join: Literal["inner", "outer"] = typer.Option(
        default="inner",
        help="How to align values when concatenating.",
    ),
    merge: Literal["same", "unique", "first", "only"] | None = typer.Option(
        default=None,
        help="Merge strategy for non-axis-aligned elements.",
    ),
    uns_merge: Literal["same", "unique", "first", "only"] | None = typer.Option(
        default=None,
        help="Merge strategy for .uns elements.",
    ),
    label: str | None = typer.Option(
        default=None,
        help="Column name in .obs/.var to store batch labels.",
    ),
    keys: list[str] = typer.Option(
        default=[],
        help="Optional keys for each object. Repeat --keys for each input file.",
    ),
    index_unique: str | None = typer.Option(
        default=None,
        help="Delimiter used to make index unique using keys.",
    ),
    fill_value: str | None = typer.Option(
        default=None,
        help="Fill value used for introduced indices with join='outer'.",
    ),
    pairwise: bool = typer.Option(
        default=False,
        help="Whether pairwise elements are included.",
    ),
):
    if not in_files:
        raise typer.BadParameter("At least one input .h5ad file is required")

    if keys and len(keys) != len(in_files):
        raise typer.BadParameter(
            f"Number of --keys values ({len(keys)}) must match number of input files ({len(in_files)})"
        )

    typer.echo(f"Concatenating {len(in_files)} files")


    ad.experimental.concat_on_disk(
        in_files=[str(p) for p in in_files],
        out_file=str(output),
        max_loaded_elems=max_loaded_elems,
        axis=axis,
        join=join,
        merge=merge,
        uns_merge=uns_merge,
        label=label,
        keys=keys or None,
        index_unique=index_unique,
        fill_value=fill_value,
        pairwise=pairwise,
    )

    typer.echo(f"Wrote concatenated AnnData to {output}")


if __name__ == "__main__":
    app()
