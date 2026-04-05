def fileSizeInGB(file) {
    def sizeGb = file.size() / (1024d * 1024d * 1024d)
    def sizeGbRoundedUp = Math.ceil(sizeGb) as long
    return "${sizeGbRoundedUp}.GB" as MemoryUnit
}

process SNAPATAC2_SPECTRAL {
    tag "$meta.id"
    container 'quay.io/cellgeni/snapatac2:2.9.0'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_spectral.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('snapatac2'), eval("python -c 'import snapatac2; print(snapatac2.__version__)'"), topic: versions, emit: versions_snapatac2

    script:
    def args = task.ext.args ?: '--n-comps 30 --features selected --random-state 0 --sample-method random --chunk-size 5000 --distance-metric cosine --weighted-by-sd --num-threads 32'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    spectral.py \\
        ${h5ad} \\
        ${prefix}_spectral.h5ad \\
        ${args}
    """
}
