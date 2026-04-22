#!/usr/bin/env python3
"""KNN neighbor graph computation that avoids touching adata.X on the GPU.

The previous version transferred the full sparse X matrix to the GPU. With a
matrix of >2^31 nonzeros this overflows int32 nnz inside cuSparse / cupyx and
silently truncates the matrix when round-tripping back to CPU, corrupting the
output h5ad.

Fix: build a slim AnnData containing only obsm[use_rep], run rapids-singlecell
neighbors on that, then merge the resulting .uns / .obsp entries into a copy
of the input file via shutil.copy + h5py.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
import cupy as cp
import h5py
import numpy as np
import pandas as pd
import anndata as ad
import rapids_singlecell as rsc

# anndata >= 0.10 moved write_elem to anndata.io; older versions had it under experimental.
try:
    from anndata.io import write_elem
except ImportError:  # pragma: no cover
    from anndata.experimental import write_elem

try:
    import hdf5plugin  # noqa: F401  -- needed for zstd-compressed snapatac2 h5ad files
except ImportError:
    pass

rmm.reinitialize(
    managed_memory=True,
    pool_allocator=False,
    devices=0,
)
cp.cuda.set_allocator(rmm_cupy_allocator)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute KNN graph with rapids-singlecell on a low-dim representation, "
                    "without touching the (potentially huge) .X matrix."
    )
    parser.add_argument("filepath", type=Path, help="Path to the input h5ad file")
    parser.add_argument("output", type=Path, help="Path to the output h5ad file")
    parser.add_argument("--n-neighbors", type=int, default=15, help="Size of local neighborhood.")
    parser.add_argument(
        "--n-pcs", type=int, default=None,
        help="Number of components from use_rep to use. None = all.",
    )
    parser.add_argument(
        "--use-rep", type=str, required=True,
        help="Representation key in obsm to use (e.g. X_spectral). REQUIRED to keep .X off the GPU.",
    )
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument(
        "--algorithm", type=str, default="brute",
        choices=["brute", "cagra", "ivfflat", "ivfpq", "nn_descent",
                 "all_neighbors", "mg_ivfflat", "mg_ivfpq"],
    )
    parser.add_argument("--metric", type=str, default="euclidean")
    parser.add_argument(
        "--method", type=str, default="umap",
        choices=["umap", "gauss", "jaccard"],
    )
    parser.add_argument(
        "--key-added", type=str, default=None,
        help="Optional key under .uns where neighbors metadata is stored.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Read input in backed mode so X never enters CPU memory or GPU memory.
    adata_in = ad.read_h5ad(args.filepath, backed="r")
    print(f"Input h5ad file:\n{adata_in}")

    if args.use_rep not in adata_in.obsm:
        raise SystemExit(f"obsm[{args.use_rep!r}] not found in input")

    rep = np.asarray(adata_in.obsm[args.use_rep])
    n_obs = adata_in.n_obs
    print(f"Using rep {args.use_rep!r} with shape {rep.shape}, dtype {rep.dtype}")

    if getattr(adata_in, "isbacked", False):
        adata_in.file.close()
    del adata_in

    # Build a slim AnnData: just obs (correct n_obs) and the chosen representation.
    # No X, no layers, no obsm/obsp/uns -> nothing big goes to the GPU.
    slim = ad.AnnData(
        obs=pd.DataFrame(index=pd.Index([str(i) for i in range(n_obs)], name="obs_id")),
        obsm={args.use_rep: rep},
    )

    # Transfer slim to GPU and compute neighbors.
    rsc.get.anndata_to_GPU(slim)
    rsc.pp.neighbors(
        slim,
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
    rsc.get.anndata_to_CPU(slim)

    print(f"Neighbors computed. New uns keys: {list(slim.uns.keys())}, "
          f"new obsp keys: {list(slim.obsp.keys())}")

    # Copy input -> output to preserve X and everything else byte-for-byte,
    # then patch in the new uns / obsp entries.
    print(f"Copying input to output: {args.filepath} -> {args.output}")
    shutil.copyfile(args.filepath, args.output)

    print("Writing neighbors entries into output via h5py")
    with h5py.File(args.output, "a") as f:
        if "obsp" not in f:
            f.create_group("obsp")
        if "uns" not in f:
            f.create_group("uns")

        for k, v in slim.obsp.items():
            if k in f["obsp"]:
                del f["obsp"][k]
            write_elem(f["obsp"], k, v)
            print(f"  wrote obsp/{k}")

        for k, v in slim.uns.items():
            if k in f["uns"]:
                del f["uns"][k]
            write_elem(f["uns"], k, v)
            print(f"  wrote uns/{k}")

    print("Done.")


if __name__ == "__main__":
    main()