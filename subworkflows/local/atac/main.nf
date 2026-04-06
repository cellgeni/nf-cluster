include { SNAPATAC2_FRAGMENTS2H5AD } from '../../../modules/local/snapatac2/fragments2h5ad'
include { AMULET } from '../../../modules/local/amulet'
include { NORMALIZECELLRANGERMETRICS } from '../../../modules/local/normalizecellrangermetrics'
include { ANNDATA_ATTACHCSV as ATTACH_CRMETRICS_AMULET } from '../../../modules/local/anndata/attachcsv'
include { ANNDATA_CONCATONDISK } from '../../../modules/local/anndata/concatondisk'
include { SNAPATAC2_METRICS } from '../../../modules/local/snapatac2/metrics'
include { SNAPATAC2_QUALITYCONTROL } from '../../../modules/local/snapatac2/qualitycontrol'
include { SNAPATAC2_ADDTILES } from '../../../modules/local/snapatac2/addtiles'
include { SNAPATAC2_SELECTFEATURES } from '../../../modules/local/snapatac2/selectfeatures'
include { SNAPATAC2_SELECTFEATURES as SNAPATAC2_SELECTFEATURES_CONCAT } from '../../../modules/local/snapatac2/selectfeatures'
include { SNAPATAC2_SCRUBLET } from '../../../modules/local/snapatac2/scrublet'
include { SNAPATAC2_SPECTRAL } from '../../../modules/local/snapatac2/spectral'
include { RAPIDS_SINGLECELL_NEIGHBORS } from '../../../modules/local/rapids_singlecell/neighbors'
include { RAPIDS_SINGLECELL_LEIDEN } from '../../../modules/local/rapids_singlecell/leiden'
include { RAPIDS_SINGLECELL_UMAP } from '../../../modules/local/rapids_singlecell/umap'

workflow ATAC {

    take:
    sample_table
    autosomes
    amulet_blacklist
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
    AMULET(amulet_input, autosomes, amulet_blacklist)

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

    // STEP5: Filter cells based on QC metrics
    SNAPATAC2_QUALITYCONTROL( SNAPATAC2_METRICS.out.h5ad )

    // STEP6: Add tile matrix
    SNAPATAC2_ADDTILES( SNAPATAC2_QUALITYCONTROL.out.h5ad )

    // STEP7: Select accessible features
    SNAPATAC2_SELECTFEATURES( SNAPATAC2_ADDTILES.out.h5ad, blacklist )

    // STEP8: Compute and filter scrublet doublets
    SNAPATAC2_SCRUBLET( SNAPATAC2_SELECTFEATURES.out.h5ad )

    // STEP9: Concatenate scrublet-processed AnnData objects on disk
    scrublet_h5ads = SNAPATAC2_SCRUBLET.out.h5ad
        .collect(flat: false)
        .map { list ->
            def ids  = list.collect { it -> it[0].id }  // Get meta from the first tuple in the list
            def h5ads = list.collect { it -> it[1] }     // Collect all h5ad paths from the list
            tuple( [id: ids], h5ads )
        }

    ANNDATA_CONCATONDISK( scrublet_h5ads )

    // STEP10: Select features on the concatenated object
    SNAPATAC2_SELECTFEATURES_CONCAT( ANNDATA_CONCATONDISK.out.h5ad, blacklist )

    // STEP11: Calculate Spectral Embedding
    SNAPATAC2_SPECTRAL( SNAPATAC2_SELECTFEATURES_CONCAT.out.h5ad )

    // STEP12: Compute KNN graph with rapids-singlecell
    RAPIDS_SINGLECELL_NEIGHBORS( SNAPATAC2_SPECTRAL.out.h5ad )

    // STEP13: Run Leiden clustering with rapids-singlecell
    RAPIDS_SINGLECELL_LEIDEN( RAPIDS_SINGLECELL_NEIGHBORS.out.h5ad )

    // STEP14: Compute UMAP embedding with rapids-singlecell
    RAPIDS_SINGLECELL_UMAP( RAPIDS_SINGLECELL_LEIDEN.out.h5ad )

    // emit:
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    // bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    // csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]
}
