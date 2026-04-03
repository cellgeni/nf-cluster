// IMPORT SUBWORKFLOWS
include { ATAC } from './subworkflows/local/atac/main'

// HELP MESSAGE
def helpMessage() {
    log.info"""
    ===========================
    Cluster ATAC-seq data
    ===========================
    Usage:
    nextflow run main.nf --sample_table sample_table.tsv
    ========================
    """.stripIndent()
}

// WORKFLOW
workflow {
    // Validate input arguments
    if (params.help) {
        helpMessage()
        System.exit(0)
    }

    // Check required arguments for peak calling
    if ( ! params.sample_table ) {
        error("Please provide --sample_table")
    }

    
    // Load files
    sample_table     = params.sample_table ? channel.value( file( params.sample_table, checkIfExists: true ) ): channel.empty()
    autosomes        = params.autosomes ? channel.value( file( params.autosomes, checkIfExists: true ) ): channel.empty()
    blacklist        = params.blacklist ? channel.value( file( params.blacklist, checkIfExists: true ) ): channel.empty()


    // Run PyCistopic pipeline
    ATAC(
        sample_table,
        autosomes,
        blacklist
    )
}