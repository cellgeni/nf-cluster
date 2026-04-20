def fileSizeInGB(file) {
    def sizeGb = file.size() / (1024d * 1024d * 1024d)
    def sizeGbRoundedUp = Math.ceil(sizeGb) as long
    return "${sizeGbRoundedUp}.GB" as MemoryUnit
}

process RAPIDS_SINGLECELL_NEIGHBORS {
    tag "KNN $meta.id"
    container 'quay.io/cellgeni/rapids-singlecell:0.14.1'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_neighbors.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('rapids-singlecell'), eval("python -c 'import rapids_singlecell; print(rapids_singlecell.__version__)'"), topic: versions, emit: versions_rapids_singlecell

    script:
    def args = task.ext.args ?: '--n-neighbors 15 --use-rep X_spectral --algorithm brute --metric euclidean --method umap --random-state 0'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    neighbors.py \\
        ${h5ad} \\
        ${prefix}_neighbors.h5ad \\
        ${args}
    """
}
