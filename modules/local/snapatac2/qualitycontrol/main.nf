def fileSizeInGB(file) {
    def sizeGb = file.size() / (1024d * 1024d * 1024d)
    def sizeGbRoundedUp = Math.ceil(sizeGb) as long
    return "${sizeGbRoundedUp}.GB" as MemoryUnit
}

process SNAPATAC2_QUALITYCONTROL {
    tag "$meta.id"
    container 'quay.io/cellgeni/snapatac2:2.9.0'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*qc.h5ad"), emit: h5ad
    tuple val(meta), path("*.png"), emit: png
    tuple val("${task.process}"), val('snapatac2'), eval("python -c 'import snapatac2; print(snapatac2.__version__)'"), topic: versions, emit: versions_snapatac2

    script:
    def args = task.ext.args ?: '--amulet-qthresh 0.01 --min-fragments 100 --max-mito-frac 0.3 --mad-thresh 3.5 --mad-key log1p_n_fragments_cr --mad-key frac_mito_reads --mad-key tsse --filter-cells'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    qualitycontrol.py ${h5ad} ${prefix} ${args}
    """
}
