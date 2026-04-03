process AMULET {
    tag "$meta.id"
    container 'quay.io/cellgeni/amulet:1.1'

    input:
    tuple val(meta), path(cellranger)
    path autosomes
    path blacklist


    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('amulet'), eval("echo v1.1"), topic: versions, emit: versions_amulet

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Check if per_barcode_metrics.csv exists in cellranger output and create singlecell.csv if it does
    if [ -f ${cellranger}/per_barcode_metrics.csv ]; then
        duckdb -c "
        COPY (
            SELECT
                barcode,
                is_cell AS is__cell_barcode
            FROM '${cellranger}/per_barcode_metrics.csv'
        ) TO 'singlecell.csv' WITH (FORMAT CSV, HEADER TRUE);
        "
    fi

    # Run AMULET
    AMULET.sh \
        ${cellranger}/*fragments.tsv.gz \
        singlecell.csv \
        ${autosomes} \
        ${blacklist} \
        .
    """
}
