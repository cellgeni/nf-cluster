# cellgeni/nf-cluster

[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

`cellgeni/nf-cluster` is a Nextflow pipeline for single-cell ATAC processing from Cell Ranger ARC output directories to integrated clustering and visualization outputs.

The current ATAC workflow performs:

- AMULET doublet calling
- metadata attachment and QC filtering
- tile matrix generation and feature selection
- Scrublet doublet scoring
- on-disk concatenation across samples
- spectral embedding
- RAPIDS neighbors, Leiden clustering, and UMAP
- Scanpy embedding plots colored by Leiden and selected metadata columns

## Usage
Prepare a sample sheet with the following columns:

```csv
sample,path
SAMPLE_A,/path/to/cellranger_arc_count_output_A
SAMPLE_B,/path/to/cellranger_arc_count_output_B
```

Each path should point to a Cell Ranger ARC output directory containing fragments files (for example, fragments.tsv.gz).

Run the ATAC workflow:

```bash
nextflow run cellgeni/nf-cluster \
   --input examples/samples.csv \
   --atac.genome hg38 \
   --outdir results
```

You can also provide parameters through a YAML/JSON params file:

```bash
nextflow run cellgeni/nf-cluster \
   -params-file params.yml
```

## Key Parameters

- Required:
   - `--input`: CSV sample sheet with sample,path
   - `--atac.genome`: genome label used by the ATAC workflow
- Common:
   - `--random_state`
   - `--outdir`
- RAPIDS neighbors:
   - `--neighbors.n_neighbors`
   - `--neighbors.algorithm`
   - `--neighbors.metric`
   - `--neighbors.method`
- RAPIDS Leiden:
   - `--leiden.resolution`
   - `--leiden.theta`
   - `--leiden.n_iterations`
   - `--leiden.key_added`
- RAPIDS UMAP:
   - `--umap.min_dist`
   - `--umap.spread`
   - `--umap.n_components`
   - `--umap.init_pos`
- Embedding plotting:
   - `--embeddingplot.basis` (default: X_umap)
   - `--embeddingplot.color` (list, default includes leiden)
   - `--embeddingplot.legend_loc`

## Output Overview

The pipeline writes linked outputs under the chosen outdir, including:

- h5ad/filtered: QC-filtered AnnData files
- h5ad/neighbors: AnnData with neighbor graph
- h5ad/leiden: AnnData with Leiden labels in obs
- h5ad/umap: AnnData with UMAP coordinates in obsm
- plots/embedding: PNG embedding plots (for example, UMAP colored by leiden)
- amulet: AMULET outputs
- reports: Nextflow execution report, timeline, trace, and DAG
