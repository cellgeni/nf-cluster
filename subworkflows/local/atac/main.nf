include { SNAPATAC2_FRAGMENTS2H5AD } from '../../../modules/local/snapatac2/fragments2h5ad'
include { AMULET } from '../../../modules/local/amulet'

workflow ATAC {

    take:
    sample_table
    autosomes
    blacklist

    main:
    // Check if genome parameter is provided
    if ( ! params.atac.genome ) {
        error("Please provide --atac.genome when using the ATAC pipeline")
    }

    // Read sample table
    cellranger = sample_table
        .splitCsv(header: true)
        .map {
            row -> tuple( [ id: row.sample ], file( row.path, checkIfExists: true ) )
        }

    // STEP1: Run AMULET on each sample
    AMULET(cellranger, autosomes, blacklist)

    // STEP2: Create H5AD files from fragments.tsv.gz
    SNAPATAC2_FRAGMENTS2H5AD( cellranger )

    // emit:
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    // bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    // csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]
}
