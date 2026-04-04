def fileSizeInGB(file) {
	def sizeGb = file.size() / (1024d * 1024d * 1024d)
	def sizeGbRoundedUp = Math.ceil(sizeGb) as long
	return "${sizeGbRoundedUp}.GB" as MemoryUnit
}

process SNAPATAC2_ADDTILES {
	tag "$meta.id"
	container 'quay.io/cellgeni/snapatac2:2.9.0'

	input:
	tuple val(meta), path(h5ad)

	output:
	tuple val(meta), path("*_tiles.h5ad"), emit: h5ad
	tuple val("${task.process}"), val('snapatac2'), eval("python -c 'import snapatac2; print(snapatac2.__version__)'"), topic: versions, emit: versions_snapatac2

	script:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${meta.id}"
	"""
	add_tiles.py \\
		${h5ad} \\
		${prefix}_tiles.h5ad \\
		${args}
	"""
}