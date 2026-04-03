process NORMALIZECELLRANGERMETRICS {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/cellgeni/duckdb:1.5.0'

    input:
    tuple val(meta), path(cellranger)

    output:
    tuple val(meta), path("*_barcodemetrics.csv"), emit: csv
    tuple val("${task.process}"), val('duckdb'), eval("duckdb --version"), topic: versions, emit: versions_normalizecellrangermetrics

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ -f ${cellranger}/per_barcode_metrics.csv ]; then
        duckdb -c "
        COPY (
            SELECT
                barcode,
                is_cell AS is__cell_barcode,
                excluded_reason,
                atac_raw_reads AS n_raw_reads,
                atac_dup_reads AS n_dup_reads,
                atac_chimeric_reads AS n_chimeric_reads,
                atac_mitochondrial_reads AS n_mitochondrial_reads,
                atac_fragments AS n_fragments,
                atac_TSS_fragments AS n_tss_fragments
            FROM '${cellranger}/per_barcode_metrics.csv'
        ) TO '${prefix}_barcodemetrics.csv' WITH (FORMAT CSV, HEADER TRUE);
        "
    elif [ -f ${cellranger}/singlecell.csv ]; then
        duckdb -c "
        COPY (
            SELECT
                barcode,
                is__cell_barcode,
                total AS n_raw_reads,
                duplicates AS n_dup_reads,
                chimeric AS n_chimeric_reads,
                mitochondrial AS n_mitochondrial_reads,
                passed_filters AS n_passed_filters,
                TSS_fragments AS n_tss_fragments,
                DNase_sensitive_region_fragments AS n_dnase_sensitive_region_fragments,
                enhancer_region_fragments AS n_enhancer_region_fragments,
                promoter_region_fragments AS n_promoter_region_fragments
            FROM '${cellranger}/singlecell.csv'
        ) TO '${prefix}_barcodemetrics.csv' WITH (FORMAT CSV, HEADER TRUE);
        "
    else
        echo "No per_barcode_metrics.csv or singlecell.csv found in ${cellranger}"
        exit 1
    fi
    """
}
