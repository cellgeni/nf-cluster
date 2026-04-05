def fileSizeInGB(file) {
    def sizeGb = file.size() / (1024d * 1024d * 1024d)
    def sizeGbRoundedUp = Math.ceil(sizeGb) as long
    return "${sizeGbRoundedUp}.GB" as MemoryUnit
}

process SNAPATAC2_SCRUBLET {
    tag "$meta.id"
    container 'quay.io/cellgeni/snapatac2:2.9.0'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_scrublet.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('snapatac2'), eval("python -c 'import snapatac2; print(snapatac2.__version__)'"), topic: versions, emit: versions_snapatac2

    script:
    def args = task.ext.args ?: '--features selected --n-comps 15 --sim-doublet-ratio 2.0 --expected-doublet-rate 0.1 --random-state 0 --probability-threshold 0.5 --no-use-approx-neighbors --n-jobs 8 --verbose'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    scrublet.py \\
        ${h5ad} \\
        ${prefix}_scrublet.h5ad \\
        ${args}
    """
}
