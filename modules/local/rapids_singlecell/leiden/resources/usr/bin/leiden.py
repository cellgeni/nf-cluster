#!/usr/bin/env python3
import argparse
from pathlib import Path

import anndata as ad
import cupy as cp
import numpy as np
import rapids_singlecell as rsc
import rmm
import scipy.sparse as sp
from rmm.allocators.cupy import rmm_cupy_allocator

rmm.reinitialize(
    managed_memory=False,  # Allows oversubscription
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
        description="Run rapids-singlecell Leiden clustering and store results in an AnnData file."
    )
    parser.add_argument("filepath", type=Path, help="Path to the input h5ad file")
    parser.add_argument("output", type=Path, help="Path to the output h5ad file")
    parser.add_argument(
        "--resolution",
        type=float,
        default=1.0,
        help="Leiden resolution parameter (higher values produce more clusters).",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=0,
        help="Random state for Leiden clustering.",
    )
    parser.add_argument(
        "--theta",
        type=float,
        default=1.0,
        help="Theta parameter for the Leiden refinement phase.",
    )
    parser.add_argument(
        "--key-added",
        type=str,
        default="leiden",
        help="obs key to store Leiden labels.",
    )
    parser.add_argument(
        "--n-iterations",
        type=int,
        default=100,
        help="Maximum number of Leiden iterations.",
    )
    parser.add_argument(
        "--use-weights",
        type=str_to_bool,
        default=True,
        help="Whether to use graph edge weights.",
    )
    parser.add_argument(
        "--neighbors-key",
        type=str,
        default=None,
        help="Neighbor metadata key in adata.uns.",
    )
    parser.add_argument(
        "--obsp",
        type=str,
        default=None,
        help="Adjacency key in adata.obsp.",
    )
    parser.add_argument(
        "--dtype",
        type=str,
        default="float32",
        choices=["float32", "float64"],
        help="Adjacency dtype used by Leiden.",
    )
    parser.add_argument(
        "--use-dask",
        type=str_to_bool,
        default=False,
        help="Use Dask multi-GPU Leiden mode.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    adata = ad.read_h5ad(args.filepath)
    print(f"Input h5ad file:\n{adata}")

    original_x_dtype = None
    cast_mode = None

    # cupyx sparse only supports '?fdFD' dtypes; cast unsupported sparse dtypes.
    if sp.issparse(adata.X) and adata.X.dtype.char not in "?fdFD":
        original_x_dtype = adata.X.dtype

        # Preserve boolean-like integer matrices by using bool on GPU and restoring
        # the original integer dtype when moving data back to CPU.
        if np.issubdtype(original_x_dtype, np.integer):
            bool_data = adata.X.data.astype(np.bool_)
            roundtrip_data = bool_data.astype(original_x_dtype, copy=False)

            if np.array_equal(roundtrip_data, adata.X.data):
                print(
                    f"Casting adata.X from {original_x_dtype} to bool for GPU compatibility"
                )
                adata.X = adata.X.astype(np.bool_)
                cast_mode = "bool"
            else:
                print(
                    f"Casting adata.X from {original_x_dtype} to float32 for GPU compatibility"
                )
                adata.X = adata.X.astype(np.float32)
                cast_mode = "float32"
        else:
            print(
                f"Casting adata.X from {original_x_dtype} to float32 for GPU compatibility"
            )
            adata.X = adata.X.astype(np.float32)
            cast_mode = "float32"

    if args.neighbors_key is not None and args.obsp is not None:
        raise ValueError("--neighbors-key and --obsp are mutually exclusive")

    # Transfer data to GPU
    rsc.get.anndata_to_GPU(adata)

    rsc.tl.leiden(
        adata,
        resolution=args.resolution,
        random_state=args.random_state,
        theta=args.theta,
        key_added=args.key_added,
        n_iterations=args.n_iterations,
        use_weights=args.use_weights,
        neighbors_key=args.neighbors_key,
        obsp=args.obsp,
        dtype=np.dtype(args.dtype),
        use_dask=args.use_dask,
        copy=False,
    )

    # Transfer results back to CPU
    rsc.get.anndata_to_CPU(adata)

    if cast_mode == "bool" and original_x_dtype is not None:
        print(f"Restoring adata.X dtype from bool to {original_x_dtype}")
        adata.X = adata.X.astype(original_x_dtype)

    print(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()
