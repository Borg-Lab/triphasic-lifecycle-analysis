# Dependencies and Upstream Provenance

This repository starts from processed RNA-seq outputs rather than raw FASTQ processing. The retained scripts reproduce the manuscript figures and supplemental tables from gene-level count/TPM matrices, annotation files, and selected downstream variant/ASE inputs.

## Upstream RNA-seq Processing

Primary expression inputs were generated before this repository by `nf-core/rnaseq` with the STAR/Salmon workflow. The upstream run metadata should be provided with the input data under:

```text
data/upstream_rnaseq/pipeline_info/
```

The recorded workflow versions in `pipeline_info/software_versions.yml` are:

```text
Nextflow: 23.10.0
nf-core/rnaseq: 3.12.0
```

The source project recorded a resumed `nf-core/rnaseq` execution using Singularity and SGE; site-specific work/output paths have been omitted from this standalone repository.

The manuscript analysis uses these nf-core/rnaseq-derived files by default, as configured in `config/paths.yml`:

```text
data/upstream_rnaseq/star_salmon/salmon.merged.gene_counts.tsv
data/upstream_rnaseq/star_salmon/salmon.merged.gene_tpm.tsv
data/upstream_rnaseq/multiqc/star_salmon/multiqc_data/multiqc_data.json
```

The nf-core/rnaseq run used STAR/Salmon and reported, among others, Salmon 1.10.1, STAR 2.7.9a/2.7.10a, FastQC 0.11.9, Trim Galore 0.6.7, Samtools 1.16.1/1.17, StringTie 2.2.1, and MultiQC-compatible pipeline reports.

## Main Manuscript Analysis Runtime

Run R scripts through the wrapper:

```bash
scripts/run_rscript.sh scripts/01_validate_inputs.R
```

By default the wrapper uses `Rscript` from `PATH`. To run inside Singularity, set these environment variables before running scripts:

```bash
export TRIPHASIC_R_IMAGE=/path/to/r-bioconductor.sif
export TRIPHASIC_R_LIBS=/path/to/R/library   # optional
export SINGULARITY_BIND=/path/to/repo,/path/to/input_data
```

The retained R scripts require R plus these CRAN/Bioconductor packages or their dependencies:

```text
AnnotationDbi, BiocGenerics, BiocManager, ComplexHeatmap, DESeq2, GO.db, GenomeInfoDb, GenomicRanges, IRanges, S4Vectors, SummarizedExperiment, WGCNA, ape, circlize, data.table, dendextend, dplyr, edgeR, ggbio, ggplot2, ggrastr, ggridges, gridExtra, myTAI, patchwork, pheatmap, purrr, readr, rlang, rtracklayer, scales, stringr, tibble, tidyr, tidyverse, yaml
```

Base/recommended R packages used by namespace include `base`, `grDevices`, `grid`, `stats`, `tools`, and `utils`.

## Allele-specific Expression Inputs and Tools

Figure 2 panel I includes an ASE preprocessing shell workflow followed by an R summarization script:

```bash
scripts/48_ase_snp_count_pipeline.sh
scripts/49_ase_expression_analysis.R
```

The ASE workflow starts from a joint VCF generated outside this repository by `nf-core/sarek`, as configured in `config/ase.yml`:

```text
data/upstream_sarek/joint_germline.vcf.gz
```

The fixed parental SNP set is built from the configured male and female VCF samples (`P1_Bm3235M_thl`, `P1_Bm3235F_thl`), then used for SNPsplit-based ASE counting. The shell workflow is written for the project cluster/SGE environment and requires:

```text
GridEngine/SGE, Singularity, conda shell initialization, bcftools, bedtools, samtools, bgzip/tabix, STAR 2.7.11b, SNPsplit 0.6.0, Trim Galore container, htseq-count container
```

The ASE script uses tool names from `PATH` by default and can be configured with `STAR_BIN`, `SNPSPLIT_BIN`, `TRIMGALORE_IMG`, `HTSEQ_IMG`, `BIND_PATHS`, and `CONDA_SH` environment variables.

## Local Configuration Files

The manuscript scripts are configured from:

```text
config/paths.yml
config/analysis.yml
config/design.tsv
config/contrasts.tsv
config/ase.yml
config/ase_samples.tsv
config/ase_sample_metadata.tsv
```

`config/paths.yml` uses repository-relative placeholder paths and writes outputs under `runs/triphasic_manuscript/results` by default. For a rerun, place local copies of the count/TPM matrices, annotation files, gene-age table, and ASE inputs at those paths or edit the config files.
