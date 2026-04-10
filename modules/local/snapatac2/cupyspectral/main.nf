def fileSizeInGB(file) {
    def sizeGb = file.size() / (1024d * 1024d * 1024d)
    sizeGb = sizeGb + Math.max(0d, sizeGb - 80d)
    def sizeGbRoundedUp = Math.ceil(sizeGb) as long
    return "${sizeGbRoundedUp}.GB" as MemoryUnit
}

process SNAPATAC2_CUPYSPECTRAL {
    tag "$meta.id"
    container 'quay.io/cellgeni/snapatac2:2.9.0'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_spectral.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('cupy'), eval("python -c 'import cupy; print(cupy.__version__)'"), topic: versions, emit: versions_cupy

    script:
    def args = task.ext.args ?: '--n-comps 30 --features selected --random-state 0 --weighted-by-sd'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    spectral.py \\
        ${h5ad} \\
        ${prefix}_spectral.h5ad \\
        ${args}
    """
}
