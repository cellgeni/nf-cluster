process SCANPY_EMBEDDING_PLOT {
    tag "$meta.id"
    container 'quay.io/cellgeni/rapids-singlecell:0.14.1'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val("${task.process}"), val('scanpy'), eval("python -c 'import scanpy as sc; print(sc.__version__)'"), topic: versions, emit: versions_scanpy

    script:
    def args = task.ext.args ?: '--basis umap --color leiden --legend-loc "right margin"'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    embeddingplot.py \\
        ${h5ad} \\
        ${prefix} \\
        ${args}
    """
}
