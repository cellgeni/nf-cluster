process RAPIDS_SINGLECELL_UMAP {
    tag "$meta.id"
    container 'quay.io/cellgeni/rapids-singlecell:0.14.1'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_umap.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('rapids-singlecell'), eval("python -c 'import rapids_singlecell; print(rapids_singlecell.__version__)'"), topic: versions, emit: versions_rapids_singlecell

    script:
    def args = task.ext.args ?: '--min-dist 0.5 --spread 1.0 --n-components 2 --alpha 1.0 --negative-sample-rate 5 --init-pos spectral --random-state 4'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    umap.py \\
        ${h5ad} \\
        ${prefix}_umap.h5ad \\
        ${args}
    """
}
