process RAPIDS_SINGLECELL_LEIDEN {
    tag "$meta.id"
    container 'quay.io/cellgeni/rapids-singlecell:0.14.1'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_leiden.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('rapids-singlecell'), eval("python -c 'import rapids_singlecell; print(rapids_singlecell.__version__)'"), topic: versions, emit: versions_rapids_singlecell

    script:
    def args = task.ext.args ?: '--resolution 1.0 --random-state 0 --theta 1.0 --key-added leiden --n-iterations 100 --use-weights true --dtype float32 --use-dask false'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    leiden.py \\
        ${h5ad} \\
        ${prefix}_leiden.h5ad \\
        ${args}
    """
}
