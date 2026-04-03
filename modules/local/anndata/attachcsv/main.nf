process ANNDATA_ATTACHCSV {
    tag "$meta.id"
    container 'quay.io/cellgeni/h5ad-cli:0.3.2'
    input:
    tuple val(meta), path(h5ad), path(metadata, stageAs: "metadata/*")

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('h5ad-cli'), eval("echo 'v0.3.2'"), topic: versions, emit: versions_anndata

    when:
    task.ext.when == null || task.ext.when

    script:
    def barcode_column = task.ext.barcode_column ?: 'barcode'
    def entry = task.ext.entry ?: 'obs'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metadata_files = metadata.collect { f -> f.toString() }.join(' ')
    """
    # Get original metadata
    h5ad export dataframe \
        ${h5ad} \
        ${entry} \
        --output metadata.csv
    
    # Merge with new metadata
    for file in ${metadata_files}; do
        duckdb -c "
            COPY (
                SELECT
                    original.* EXCLUDE (_index),
                    new.*
                FROM 'metadata.csv' AS original
                LEFT JOIN '\$file' AS new 
                    ON original._index = new.${barcode_column}
            ) TO 'metadata_merged.csv' WITH (FORMAT CSV, HEADER TRUE);
        "
        mv metadata_merged.csv metadata.csv
    done

    
    h5ad import dataframe \
        ${h5ad} \
        ${entry} \
        metadata.csv \
        --output ${prefix}.h5ad \
        --index-column ${barcode_column}
    """
}
