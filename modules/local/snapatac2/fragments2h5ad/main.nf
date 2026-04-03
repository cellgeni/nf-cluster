process SNAPATAC2_FRAGMENTS2H5AD {
    tag "$meta.id"
    container 'quay.io/cellgeni/snapatac2:2.9.0'

    input:
    tuple val(meta), path(cellranger)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('snapatac2'), eval("python -c 'import snapatac2; print(snapatac2.__version__)'"), topic: versions, emit: versions_snapatac2

    script:
    def args = task.ext.args ?: '--genome hg38 --min-num-fragments 100'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filtered = task.ext.filtered || true 
    """
    fragments=${cellranger}/*fragments.tsv.gz
    whitelist=${cellranger}/filtered*matrix/barcodes.tsv*

    # check if fragments file exists
    if [ ! -f \$fragments ]; then
        echo "No fragments file found in ${cellranger}"
        exit 1
    fi

    # check if whitelist file exists if filtered is true
    if [ "${filtered}" = "true" ] && [ ! -f \$whitelist ]; then
        echo "No whitelist file found in ${cellranger} but --filtered is set to true"
        exit 1
    fi

    fragments2h5ad.py \
        \$fragments \
        ${prefix}.h5ad \
        --sample-key "${meta.id}" \
        $args \
        ${filtered ? "--whitelist \$whitelist" : ""} \
        --n-jobs ${task.cpus}
    """
}
