#!/usr/bin/env python3
import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot embedding panels from an AnnData file using scanpy.pl.embedding."
    )
    parser.add_argument("filepath", type=Path, help="Path to the input h5ad file")
    parser.add_argument("output_prefix", type=Path, help="Prefix for output image files")
    parser.add_argument(
        "--basis",
        type=str,
        default="umap",
        help="Embedding basis name, e.g. umap for adata.obsm['X_umap'].",
    )
    parser.add_argument(
        "--color",
        action="append",
        default=["leiden"],
        help="Observation or gene key to color by. Repeat for multiple panels.",
    )
    parser.add_argument(
        "--legend-loc",
        type=str,
        default="right margin",
        help="Legend location passed to scanpy.pl.embedding.",
    )
    parser.add_argument(
        "--ncols",
        type=int,
        default=4,
        help="Number of panel columns when plotting multiple color keys.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="PNG output DPI.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    adata = ad.read_h5ad(args.filepath)
    print(f"Input h5ad file:\n{adata}")

    if args.basis not in adata.obsm:
        raise ValueError(f"Embedding basis '{args.basis}' not found in adata.obsm[{args.basis!r}]")

    colors = []
    seen = set()
    for key in args.color:
        if key in seen:
            continue
        seen.add(key)
        if key in adata.obs.columns or key in adata.var_names:
            colors.append(key)
        else:
            print(f"Skipping color key '{key}' because it was not found in adata.obs or adata.var_names")

    color_arg = colors if colors else None

    sc.pl.embedding(
        adata,
        basis=args.basis,
        color=color_arg,
        legend_loc=args.legend_loc,
        ncols=args.ncols,
        show=False,
    )

    output_path = Path(f"{args.output_prefix}_{args.basis.replace('X_', '')}.png")
    plt.savefig(output_path, bbox_inches="tight", dpi=args.dpi)
    plt.close("all")

    print(f"Saved embedding plot: {output_path}")


if __name__ == "__main__":
    main()
