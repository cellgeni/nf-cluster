#!/usr/bin/env python3
import argparse
from pathlib import Path

import anndata as ad
import cupy as cp
import rapids_singlecell as rsc
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

rmm.reinitialize(
    managed_memory=True,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)


def str_to_bool(value: str) -> bool:
    if isinstance(value, bool):
        return value

    value_l = value.lower()
    if value_l in {"true", "1", "yes", "y"}:
        return True
    if value_l in {"false", "0", "no", "n"}:
        return False

    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run rapids-singlecell UMAP embedding and store results in an AnnData file."
    )
    parser.add_argument("filepath", type=Path, help="Path to the input h5ad file")
    parser.add_argument("output", type=Path, help="Path to the output h5ad file")
    parser.add_argument("--min-dist", type=float, default=0.5, help="UMAP min_dist parameter.")
    parser.add_argument("--spread", type=float, default=1.0, help="UMAP spread parameter.")
    parser.add_argument(
        "--n-components",
        type=int,
        default=2,
        help="Number of embedding dimensions.",
    )
    parser.add_argument(
        "--maxiter",
        type=int,
        default=None,
        help="Maximum number of UMAP optimization iterations.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=1.0,
        help="Initial learning rate for UMAP optimization.",
    )
    parser.add_argument(
        "--negative-sample-rate",
        type=int,
        default=5,
        help="Number of negative samples per positive sample.",
    )
    parser.add_argument(
        "--init-pos",
        type=str,
        default="auto",
        help="Initialization strategy for UMAP (auto/spectral/random/paga or obsm key).",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=0,
        help="Random state for UMAP embedding.",
    )
    parser.add_argument(
        "--a",
        type=float,
        default=None,
        help="UMAP a parameter. If unset, inferred from min_dist/spread.",
    )
    parser.add_argument(
        "--b",
        type=float,
        default=None,
        help="UMAP b parameter. If unset, inferred from min_dist/spread.",
    )
    parser.add_argument(
        "--key-added",
        type=str,
        default=None,
        help="Optional key to store UMAP output in obsm/uns.",
    )
    parser.add_argument(
        "--neighbors-key",
        type=str,
        default=None,
        help="Optional neighbors key in adata.uns.",
    )
    parser.add_argument(
        "--copy",
        type=str_to_bool,
        default=False,
        help="Return a copy instead of modifying adata in place.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    adata = ad.read_h5ad(args.filepath)
    print(f"Input h5ad file:\n{adata}")

    # Transfer data to GPU
    rsc.get.anndata_to_GPU(adata)

    rsc.tl.umap(
        adata,
        min_dist=args.min_dist,
        spread=args.spread,
        n_components=args.n_components,
        maxiter=args.maxiter,
        alpha=args.alpha,
        negative_sample_rate=args.negative_sample_rate,
        init_pos=args.init_pos,
        random_state=args.random_state,
        a=args.a,
        b=args.b,
        key_added=args.key_added,
        neighbors_key=args.neighbors_key,
        copy=args.copy,
    )

    # Transfer results back to CPU
    rsc.get.anndata_to_CPU(adata)

    print(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()
