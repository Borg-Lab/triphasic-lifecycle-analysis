# Manuscript Code: Regulatory evolution and architecture of the triphasic life cycle in the red alga *Bostrychia moritziana*  

This repository contains the code used for transcriptomic analysis of the triphasic life history of the model red alga *Bostrychia moritziana*.

## Scope

This repository is the code archive for the manuscript. It contains analysis scripts, configuration files, and final figure and table outputs.

## Repository Structure

- `scripts/`: ordered analysis scripts and the complete pipeline runner
- `R/`: shared R helper functions used in the analysis scripts
- `config/`: analysis configuration, sample design, and contrasts
- `manuscript_outputs/`: final manuscript figures and supplemental tables
- `metadata/`: manifests mapping each manuscript output to its generating script
- `docs/`: dependency and data-provenance documentation

## Inputs

The paths in `config/paths.yml` are repository-relative placeholders for the required expression matrices and annotation files. Place the corresponding files at those locations or update the configuration before running the analysis. ASE analysis has additional inputs defined in `config/ase.yml`.

The required software and upstream data provenance are documented in `docs/dependencies_and_provenance.md`.

## Running the Analysis

The final figures and tables used for the manuscript are already provided under `manuscript_outputs/`; rerunning the pipeline is not required to inspect them.

To rebuild the analyses from the configured inputs, run the complete pipeline from the repository root:

```bash
scripts/run_manuscript_pipeline.sh
```

The runner executes the numbered scripts in dependency order, including the computationally intensive ASE workflow. Individual numbered scripts may be run separately after their prerequisite outputs have been generated, for example:

```bash
scripts/run_rscript.sh scripts/30_plot_tai_signature_replots.R
```

The following environment variables can be used to select alternative configuration files:

```bash
PATHS_CONFIG=config/paths.yml
ANALYSIS_CONFIG=config/analysis.yml
DESIGN_CONFIG=config/design.tsv
CONTRASTS_CONFIG=config/contrasts.tsv
```

By default, generated results are written to `runs/triphasic_manuscript/results`.

## Output Mapping

`metadata/manuscript_output_manifest.tsv` maps each manuscript panel and supplemental table to its archived output, source result, and generating script. This is the primary reference for identifying the code associated with a particular figure.

## Upstream Processing

This manuscript analysis uses processed STAR/Salmon outputs from the RNA-seq processing pipeline `nf-core/rnaseq` 3.12.0 and Nextflow 23.10.0.

## License

The code in this repository is released under the MIT License. See `LICENSE`.
