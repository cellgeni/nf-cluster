#!/usr/bin/env python3
from pathlib import Path
from typing import Annotated, Literal

import anndata as ad
import scipy.sparse as sp

import cupy as cp
import cupyx.scipy.sparse as cpx_sparse
from cupyx.scipy.sparse.linalg import LinearOperator, eigsh
import typer
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

rmm.reinitialize(
    managed_memory=True,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)


app = typer.Typer()

_ROW_SCALE_KERNEL = cp.RawKernel(
    r"""
extern "C" __global__
void row_scale(const int n_rows, const int* indptr, float* data, const float* scales) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row < n_rows) {
        const float s = scales[row];
        for (int j = indptr[row]; j < indptr[row + 1]; ++j) {
            data[j] *= s;
        }
    }
}
""",
    "row_scale",
)


def scale_rows_inplace(X_csr, scales):
    """Multiply each row by its scale, in place."""
    scales = cp.asarray(scales, dtype=cp.float32)
    n_rows = X_csr.shape[0]
    blocks = ((n_rows + 255) // 256,)
    _ROW_SCALE_KERNEL(blocks, (256,), (n_rows, X_csr.indptr, X_csr.data, scales))
    return X_csr


def row_sq_sum(X_csr):
    """Row-wise sum of squares for CSR without sparse-sparse multiply."""
    data_sq = X_csr.data * X_csr.data
    prefix = cp.empty(data_sq.size + 1, dtype=cp.float32)
    prefix[0] = 0.0
    cp.cumsum(data_sq, out=prefix[1:])
    return prefix[X_csr.indptr[1:]] - prefix[X_csr.indptr[:-1]]


def to_gpu_csr(X):
    """Convert scipy/cupyx sparse matrix to cupyx CSR float32."""
    if cpx_sparse.isspmatrix(X):
        X = X.tocsr()
        X.sum_duplicates()
        X.eliminate_zeros()
        return X.astype(cp.float32, copy=False)

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

    raise TypeError("X must be a scipy.sparse or cupyx.scipy.sparse matrix")


def idf_gpu(X_csr):
    """
    SnapATAC2-style IDF:
        log(n_cells / (1 + df))
    where df is the number of cells with a nonzero count per feature.
    """
    n_cells = X_csr.shape[0]
    n_cells, n_features = X_csr.shape
    df = cp.bincount(X_csr.indices, minlength=n_features).astype(cp.float32, copy=False)
    return cp.log(cp.float32(n_cells) / (1.0 + df))


def scale_columns_inplace(X_csr, weights):
    """
    Multiply each feature column by its weight.
    Equivalent to preprocessing by IDF before cosine similarity.
    """
    X_csr = X_csr.copy()
    X_csr.data *= weights[X_csr.indices]
    return X_csr


def row_l2_normalize(X_csr, eps=1e-12):
    """Return row-L2-normalized sparse matrix."""
    row_sq = cp.asarray(X_csr.multiply(X_csr).sum(axis=1)).ravel()
    row_norm = cp.sqrt(cp.maximum(row_sq, eps))
    inv_row = 1.0 / row_norm
    return scale_rows_inplace(X_csr, inv_row)


def snapatac2_cosine_spectral_gpu(
    X,
    n_comps=30,
    feature_weights=None,  # if None, use SnapATAC2-style IDF
    weighted_by_sd=True,
    random_state=0,
    eps=1e-8,
):
    """
    Closest CuPy/CuPyX match to SnapATAC2 cosine spectral embedding.

    Parameters
    ----------
    X : scipy.sparse or cupyx sparse matrix
        Cell-by-feature count matrix.
    features : None, boolean mask, or integer indices
        Selected features. Apply selection BEFORE IDF, matching Snap usage.
    n_comps : int
        Number of eigenvectors to request.
    feature_weights : array-like or None
        If None, uses IDF.
    weighted_by_sd : bool
        If True, keep only positive eigenvalues and multiply eigenvectors
        by sqrt(eigenvalues), matching SnapATAC2 behavior.
    random_state : int
        Seed for eigensolver start vector.

    Returns
    -------
    evals, evecs : cupy.ndarray, cupy.ndarray
    """
    X = to_gpu_csr(X)

    n_cells, n_features = X.shape
    if n_cells < 2:
        raise ValueError("Need at least 2 cells")
    k = min(int(n_comps), n_cells - 1)
    if k < 1:
        raise ValueError("n_comps must be >= 1")

    # Snap docs describe preprocessing by IDF before cosine similarity.
    if feature_weights is None:
        w = idf_gpu(X)
    else:
        w = cp.asarray(feature_weights, dtype=cp.float32)
        if w.shape[0] != n_features:
            raise ValueError(
                "feature_weights length must match number of selected features"
            )

    # 1) IDF weighting
    Xw = scale_columns_inplace(X, w)

    # 2) Row L2 normalization -> X_norm
    Xn = row_l2_normalize(Xw)

    # 3) Degree for W = X_norm X_norm^T - I
    #    degree = (X_norm @ (X_norm.T @ 1)) - 1
    ones = cp.ones(n_cells, dtype=cp.float32)
    xt1 = Xn.T @ ones
    degree = cp.asarray(Xn @ xt1).ravel() - 1.0
    degree = cp.maximum(degree, eps)

    dinv_sqrt = 1.0 / cp.sqrt(degree)
    dinv = 1.0 / degree

    # \tilde{X} = D^{-1/2} X_norm
    Xtilde = scale_rows_inplace(Xn, dinv_sqrt)

    # Matrix-free operator:
    # \tilde{W} v = \tilde{X} (\tilde{X}^T v) - D^{-1} v
    def matvec(v):
        v = cp.asarray(v, dtype=cp.float32)
        return Xtilde @ (Xtilde.T @ v) - dinv * v

    op = LinearOperator(
        shape=(n_cells, n_cells),
        matvec=matvec,
        dtype=cp.float32,
    )

    # Reproducible start vector
    rs = cp.random.RandomState(random_state)
    v0 = rs.standard_normal(n_cells).astype(cp.float32)
    v0 /= cp.linalg.norm(v0) + eps

    # Match Snap's fallback behavior more closely with LM
    evals, evecs = eigsh(op, k=k, which="LM", v0=v0)

    # Sort descending like the Snap code
    order = cp.argsort(evals)[::-1]
    evals = cp.real(evals[order])
    evecs = cp.real(evecs[:, order])

    # Match Snap's weighted_by_sd postprocessing
    if weighted_by_sd:
        keep = evals > 0
        evals = evals[keep]
        evecs = evecs[:, keep] * cp.sqrt(evals)

    return evals, evecs


@app.command()
def main(
    filepath: Annotated[Path, typer.Argument(help="Path to the input h5ad file")],
    output: Annotated[Path, typer.Argument(help="Path to the output h5ad file")],
    n_comps: int = typer.Option(
        default=30,
        help="Number of dimensions to keep.",
    ),
    features: str = typer.Option(
        default="selected",
        help="Feature mask name in adata.var, or None to use all features.",
    ),
    random_state: int = typer.Option(
        default=0,
        help="Seed of the random state generator.",
    ),
    weighted_by_sd: bool = typer.Option(
        default=True,
        help="Whether to weight eigenvectors by square root of eigenvalues.",
    ),
    feature_weights: list[float] = typer.Option(
        default=[],
        help="Optional feature weights; repeat --feature-weights to pass multiple values.",
    ),
):
    adata = ad.read_h5ad(filepath)
    typer.echo(f"Input h5ad file:\n{adata}")

    parsed_features = None if features.lower() == "none" else features

    # Subset features if a mask is provided
    X = (
        adata[:, adata.var[parsed_features]].X.copy()
        if parsed_features is not None
        else adata.X
    )

    typer.echo(f"selected_features = {int(adata.var['selected'].sum())}")
    typer.echo(f"shape = {X.shape}")
    typer.echo(f"nnz = {int(X.nnz)}")
    typer.echo(f"nnz fits int32 = {int(X.nnz) < 2_147_483_647}")

    evals, evecs = snapatac2_cosine_spectral_gpu(
        X,
        n_comps=n_comps,
        feature_weights=feature_weights if feature_weights else None,
        weighted_by_sd=weighted_by_sd,
        random_state=random_state,
    )

    adata.obsm["X_spectral"] = evecs.get()
    adata.uns["spectral_eigenvalue"] = evals.get()

    typer.echo(f"Output h5ad file:\n{adata}")
    adata.write_h5ad(output)


if __name__ == "__main__":
    app()
