def fileSizeInGB(file) {
    def sizeGb = file.size() / (1024d * 1024d * 1024d)
    def sizeGbRoundedUp = Math.ceil(sizeGb) as long
    return "${sizeGbRoundedUp}.GB" as MemoryUnit
}

process RAPIDS_SINGLECELL_SELECTFEATURES {
    tag "$meta.id"
    container 'quay.io/cellgeni/snapatac2:2.9.0'

    input:
    tuple val(meta), path(h5ad)
    path blacklist

    output:
    tuple val(meta), path("*_features.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('rapids-singlecell'), eval("python -c 'import rapids_singlecell; print(rapids_singlecell.__version__)'"), topic: versions, emit: versions_rapids_singlecell

    script:
    def args = task.ext.args ?: '--n-features 50000 --random-state 0'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_features.py \\
        ${h5ad} \\
        ${prefix}_features.h5ad \\
        --blacklist ${blacklist} \\
        ${args}
    """
}
