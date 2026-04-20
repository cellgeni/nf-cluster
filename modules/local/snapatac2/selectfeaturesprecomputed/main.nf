def fileSizeInGB(file) {
    def sizeGb = file.size() / (1024d * 1024d * 1024d)
    def sizeGbRoundedUp = Math.ceil(sizeGb) as long
    return "${sizeGbRoundedUp}.GB" as MemoryUnit
}

process SNAPATAC2_SELECTFEATURESPRECOMPUTED {
    tag "$meta.id"
    container 'quay.io/cellgeni/snapatac2:2.9.0'

    input:
    tuple val(meta), path(h5ad)
    path blacklist

    output:
    tuple val(meta), path("*_features.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('snapatac2'), eval("python -c 'import snapatac2; print(snapatac2.__version__)'"), topic: versions, emit: versions_snapatac2

    script:
    def args = task.ext.args ?: '--groupby leiden_selectfeatures --n-iter 1 --n-features 50000 --subset-selected'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_features.py \\
        ${h5ad} \\
        ${prefix}_features.h5ad \\
        --blacklist ${blacklist} \\
        ${args}
    """
}
