#!/usr/bin/env python3
import argparse
from pathlib import Path

import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
import cupy as cp
import numpy as np
import anndata as ad
import rapids_singlecell as rsc
import scipy.sparse as sp

rmm.reinitialize(
    managed_memory=False,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute KNN graph with rapids-singlecell and store results in an AnnData file."
    )
    parser.add_argument("filepath", type=Path, help="Path to the input h5ad file")
    parser.add_argument("output", type=Path, help="Path to the output h5ad file")
    parser.add_argument("--n-neighbors", type=int, default=15, help="Size of local neighborhood.")
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=None,
        help="Number of principal components to use. If None, rapids-singlecell defaults are used.",
    )
    parser.add_argument(
        "--use-rep",
        type=str,
        default=None,
        help="Representation key to use (e.g. X_spectral).",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=0,
        help="Random state for neighbor graph computation.",
    )
    parser.add_argument(
        "--algorithm",
        type=str,
        choices=[
            "brute",
            "cagra",
            "ivfflat",
            "ivfpq",
            "nn_descent",
            "all_neighbors",
            "mg_ivfflat",
            "mg_ivfpq",
        ],
        default="brute",
        help="Neighbor search algorithm.",
    )
    parser.add_argument("--metric", type=str, default="euclidean", help="Distance metric.")
    parser.add_argument(
        "--method",
        type=str,
        choices=["umap", "gauss", "jaccard"],
        default="umap",
        help="Method for computing connectivities.",
    )
    parser.add_argument(
        "--key-added",
        type=str,
        default=None,
        help="Optional key under .uns where neighbors metadata is stored.",
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

    # Transfer data to GPU
    rsc.get.anndata_to_GPU(adata)

    # Compute neighbors using rapids-singlecell
    rsc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        n_pcs=args.n_pcs,
        use_rep=args.use_rep,
        random_state=args.random_state,
        algorithm=args.algorithm,
        metric=args.metric,
        method=args.method,
        key_added=args.key_added,
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
