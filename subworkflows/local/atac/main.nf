include { SNAPATAC2_FRAGMENTS2H5AD } from '../../../modules/local/snapatac2/fragments2h5ad'
include { AMULET } from '../../../modules/local/amulet'
include { NORMALIZECELLRANGERMETRICS } from '../../../modules/local/normalizecellrangermetrics'
include { ANNDATA_ATTACHCSV as ATTACH_CRMETRICS_AMULET } from '../../../modules/local/anndata/attachcsv'
include { SNAPATAC2_METRICS } from '../../../modules/local/snapatac2/metrics'

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
    sample_table
        .splitCsv(header: true)
        .map {
            row -> tuple( [ id: row.sample ], file( row.path, checkIfExists: true ) )
        }
        .tap { cellranger }
        .map {
            meta, path -> 
            def fragments = file( "${path.toString()}/*fragments.{tsv,tsv.gz}" )[0]
            if ( ! fragments.exists() ) {
                    error("No fragments file found for sample ${meta.id} at path: ${fragments}. Please check your sample table.")
                }
            tuple( meta, fragments )
        }
        .set { fragments }

    // STEP0: Normalize CellRanger metrics
    NORMALIZECELLRANGERMETRICS( cellranger )

    // STEP1: Run AMULET on each sample
    amulet_input = fragments.join( NORMALIZECELLRANGERMETRICS.out.csv, by: [0], failOnDuplicate: true, failOnMismatch: true )
    AMULET(amulet_input, autosomes, blacklist)

    // STEP2: Create H5AD files from fragments.tsv.gz
    SNAPATAC2_FRAGMENTS2H5AD( cellranger )

    // STEP3: Attach AMULET results and normalized CellRanger metrics to H5AD files
    // Collect AMULET results and normalized CellRanger metrics for each sample
    metadata = AMULET.out.csv
        .mix( 
            NORMALIZECELLRANGERMETRICS.out.csv,
            )
        .groupTuple(by: [0], size: 2, remainder: true)
    
    // Join metadata with H5AD files
    h5ad_metadata = SNAPATAC2_FRAGMENTS2H5AD.out.h5ad
        .join( metadata, by: [0], failOnDuplicate: true, failOnMismatch: true )

    ATTACH_CRMETRICS_AMULET(h5ad_metadata)

    // STEP4: Calculate additional metrics using SNAPATAC2
    SNAPATAC2_METRICS( ATTACH_CRMETRICS_AMULET.out.h5ad )

    // emit:
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    // bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    // csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]
}
