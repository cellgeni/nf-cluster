def filesSizeInGB(file) {
	def totalSize = file.collect { it.size() }.sum()
    def sizeGb = totalSize / (1024d * 1024d * 1024d)
    def sizeGbRoundedUp = Math.ceil(sizeGb / 10 ) as long
    return "${sizeGbRoundedUp}.GB" as MemoryUnit
}

process ANNDATA_CONCATONDISK {
    tag "Concatenating ${meta.id.size()} files"
    container 'quay.io/cellgeni/snapatac2:2.9.0'

    input:
    tuple val(meta), path(h5ads)

    output:
    tuple val(new_meta), path("*.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('anndata'), eval("python -c 'import anndata; print(anndata.__version__)'"), topic: versions, emit: versions_anndata

    script:
    def args = task.ext.args ?: "--axis obs --join inner --merge same --uns-merge same --index-unique '___'"
    def keys = meta.id instanceof List ? "--keys ${meta.id.join(' --keys ')}" : ""
    new_meta = [id: "concatenated", oldid: meta.id]
    """
    concat_on_disk.py \\
        ${h5ads} \\
        ${keys} \\
        --output concatenated.h5ad ${args}
    """
}
