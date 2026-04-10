#!/usr/bin/env python3
import gc
from typing import Annotated
from pathlib import Path

import typer
import anndata as ad
import numpy as np
import scipy.sparse as sp

import snapatac2._snapatac2 as snap_internal

import cupy as cp
import cupyx.scipy.sparse as cpx_sparse
from cupyx.scipy.sparse.linalg import LinearOperator, eigsh

import rapids_singlecell as rsc
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

rmm.reinitialize(
    managed_memory=True,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)


app = typer.Typer()

def free_gpu_memory():
    # 1) Drop Python references first
    gc.collect()

    # 2) Make sure queued CUDA work is finished
    cp.cuda.Stream.null.synchronize()

    # 3) If you are using CuPy's default allocator, release cached free blocks
    cp.get_default_memory_pool().free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()


def to_gpu_csr(X):
    """Convert scipy/cupyx sparse matrix to cupyx CSR float32."""
    if cpx_sparse.isspmatrix(X):
        X = X.tocsr()
        X.sum_duplicates()
        X.eliminate_zeros()
        return X.astype(cp.float32)

    if sp.issparse(X):
        X = X.tocsr()
        X.sum_duplicates()
        X.eliminate_zeros()
        return cpx_sparse.csr_matrix(
            (
                cp.asarray(X.data, dtype=cp.float32),
                cp.asarray(X.indices, dtype=cp.int32),
                cp.asarray(X.indptr, dtype=cp.int32),
            ),
            shape=X.shape,
        )

    raise TypeError(f"X must be a scipy.sparse or cupyx.scipy.sparse matrix, got {type(X)}")


def idf_gpu(X_csr):
    """
    SnapATAC2-style IDF:
        log(n_cells / (1 + df))
    where df is the number of cells with a nonzero count per feature.
    """
    n_cells = X_csr.shape[0]
    X_csc = X_csr.tocsc(copy=True)
    df = cp.diff(X_csc.indptr).astype(cp.float32)
    return cp.log(cp.float32(n_cells) / (1.0 + df))


def scale_columns_inplace(X_csr, weights):
    """
    Multiply each feature column by its weight.
    """
    X_csr = X_csr.copy()
    X_csr.data *= weights[X_csr.indices]
    return X_csr


def row_l2_normalize(X_csr, eps=1e-12):
    """Return row-L2-normalized sparse matrix."""
    row_sq = cp.asarray(X_csr.multiply(X_csr).sum(axis=1)).ravel()
    row_norm = cp.sqrt(cp.maximum(row_sq, eps))
    inv_row = 1.0 / row_norm
    return cpx_sparse.diags(inv_row) @ X_csr


def snapatac2_cosine_spectral_gpu(
    X,
    features=None,              # boolean mask or integer indices
    n_comps=30,
    feature_weights=None,       # if None, use SnapATAC2-style IDF
    weighted_by_sd=True,
    random_state=0,
    eps=1e-8,
):
    """
    CuPy/CuPyX exact cosine spectral embedding matching SnapATAC2's visible behavior.
    """
    X = to_gpu_csr(X)

    if features is not None:
        features = cp.asarray(features)
        X = X[:, features]

    n_cells, n_features = X.shape
    if n_cells < 2:
        raise ValueError("Need at least 2 cells")
    k = min(int(n_comps), n_cells - 1, n_features - 1)
    if k < 1:
        raise ValueError("n_comps must be >= 1")

    if feature_weights is None:
        w = idf_gpu(X)
    else:
        w = cp.asarray(feature_weights, dtype=cp.float32)
        if w.shape[0] != n_features:
            raise ValueError("feature_weights length must match number of selected features")

    Xw = scale_columns_inplace(X, w)
    Xn = row_l2_normalize(Xw)

    ones = cp.ones(n_cells, dtype=cp.float32)
    xt1 = Xn.T @ ones
    degree = cp.asarray(Xn @ xt1).ravel() - 1.0
    degree = cp.maximum(degree, eps)

    dinv_sqrt = 1.0 / cp.sqrt(degree)
    dinv = 1.0 / degree

    Xtilde = cpx_sparse.diags(dinv_sqrt) @ Xn

    def matvec(v):
        v = cp.asarray(v, dtype=cp.float32)
        return Xtilde @ (Xtilde.T @ v) - dinv * v

    op = LinearOperator(
        shape=(n_cells, n_cells),
        matvec=matvec,
        dtype=cp.float32,
    )

    rs = cp.random.RandomState(random_state)
    v0 = rs.standard_normal(n_cells).astype(cp.float32)
    v0 /= cp.linalg.norm(v0) + eps

    evals, evecs = eigsh(op, k=k, which="LM", v0=v0)

    order = cp.argsort(evals)[::-1]
    evals = cp.real(evals[order])
    evecs = cp.real(evecs[:, order])

    if weighted_by_sd:
        keep = evals > 0
        evals = evals[keep]
        evecs = evecs[:, keep] * cp.sqrt(evals)

    return evals, evecs


def _find_most_accessible_features_gpu(
    feature_count,
    filter_lower_quantile,
    filter_upper_quantile,
    total_features,
):
    idx = cp.argsort(feature_count)
    sorted_counts = feature_count[idx]
    idx = idx[sorted_counts > 0]

    n = int(idx.size)
    if n == 0:
        return idx

    n_lower = int(filter_lower_quantile * n)
    n_upper = int(filter_upper_quantile * n)

    if n_lower + n_upper >= n:
        return cp.empty((0,), dtype=cp.int64)

    idx = idx[n_lower:n - n_upper]
    return idx[::-1][:total_features]


def _blacklist_mask_cpu(var_names, blacklist_path: Path | None):
    if blacklist_path is None:
        return None
    return np.asarray(snap_internal.intersect_bed(var_names, str(blacklist_path)), dtype=bool)


def select_features_gpu(
    adata,
    n_features=50000,
    blacklist: Path | None = None,
    random_state=0,
    neighbors_k=15,
    leiden_resolution=1.0,
):
    """
    GPU-first reimplementation of SnapATAC2 select_features.
    """
    
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
                typer.echo(
                        f"Casting adata.X from {original_x_dtype} to bool for GPU compatibility"
                    )
                adata.X = adata.X.astype(np.bool_)
                cast_mode = "bool"
            else:
                typer.echo(
                        f"Casting adata.X from {original_x_dtype} to float32 for GPU compatibility"
                    )
                adata.X = adata.X.astype(np.float32)
                cast_mode = "float32"
        else:
            
            typer.echo(
                    f"Casting adata.X from {original_x_dtype} to float32 for GPU compatibility"
                )
            adata.X = adata.X.astype(np.float32)
            cast_mode = "float32"

    typer.echo(f"Move AnnData to GPU for neighbors and Leiden...")
    n_obs, n_vars = adata.shape
    rsc.get.anndata_to_GPU(adata)

    
    typer.echo(f"Computing neighbors and Leiden clustering for feature selection refinement...")

    rsc.pp.neighbors(
        adata,
        use_rep="X_spectral",
        n_neighbors=neighbors_k,
        algorithm="brute",
        metric="euclidean",
        random_state=random_state,
    )

    
    typer.echo(f"Running Leiden clustering for feature selection refinement...")

    rsc.tl.leiden(
        adata,
        resolution=leiden_resolution,
        key_added="_gpu_leiden",
        random_state=random_state,
    )

    
    typer.echo(f"Aggregating counts by Leiden clusters and selecting features with highest variance...")
    agg = rsc.get.aggregate(
        adata,
        by="_gpu_leiden",
        func="sum",
        return_sparse=False,
    )

    typer.echo(f"Aggregated anndata: {agg}")
    typer.echo(f"Finding most accessible features after iteration...")

    var = cp.var(cp.log1p(agg.layers["sum"]), axis=0)
    selected_features = cp.argsort(var)[::-1][:n_features]

    blacklist_mask_cpu = _blacklist_mask_cpu(adata.var_names, blacklist)
    if blacklist_mask_cpu is not None and selected_features.size > 0:
        blacklist_mask = cp.asarray(blacklist_mask_cpu)
        selected_features = selected_features[~blacklist_mask[selected_features]]

    
    typer.echo(f"Selected features after iteration: {selected_features.size} out of {n_vars}")


    result = cp.zeros(n_vars, dtype=cp.bool_)
    result[selected_features] = True

    adata.var["selected"] = cp.asnumpy(result)

    if "_gpu_leiden" in adata.obs:
        del adata.obs["_gpu_leiden"]
    if "_gpu_leiden" in adata.uns:
        del adata.uns["_gpu_leiden"]
    if "X_spectral" in adata.obsm:
        del adata.obsm["X_spectral"]
    if "neighbors" in adata.uns:
        adata.uns.pop("neighbors", None)
    if "distances" in adata.obsp:
        adata.obsp.pop("distances", None)
    if "connectivities" in adata.obsp:
        adata.obsp.pop("connectivities", None)

    rsc.get.anndata_to_CPU(adata)

    if cast_mode == "bool" and original_x_dtype is not None:
        typer.echo(f"Restoring adata.X dtype from bool to {original_x_dtype}")
        adata.X = adata.X.astype(original_x_dtype)
    elif cast_mode == "float32" and original_x_dtype is not None:
        typer.echo(f"Restoring adata.X dtype from float32 to {original_x_dtype}")
        adata.X = adata.X.astype(original_x_dtype)


@app.command()
def main(
    filepath: Annotated[Path, typer.Argument(help="Path to the input h5ad file")],
    output: Annotated[Path, typer.Argument(help="Path to the output h5ad file")],
    n_features: int = typer.Option(
        default=50000,
        help="Number of features to select based on accessibility.",
    ),
    blacklist: Path | None = typer.Option(
        default=None,
        help="Path to a BED file containing genomic regions to exclude from feature selection.",
    ),
    random_state: int = typer.Option(
        default=0,
        help="Random seed used by GPU spectral embedding and Leiden.",
    ),
    subset_selected: bool = typer.Option(
        default=False,
        help="Whether to subset the AnnData object to only include the selected features. If False, the selected features will be added to .var['selected'] and the original feature matrix will be retained.",
    ),
):
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    select_features_gpu(
        adata,
        n_features=n_features,
        blacklist=blacklist,
        random_state=random_state,
    )

    if subset_selected:
        adata = adata[:, adata.var["selected"].to_numpy()].copy()

    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()
