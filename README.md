# Manuscript Code: Regulatory evolution and architecture of the triphasic life cycle in the red alga *Bostrychia moritziana*  

This repository contains the code used for transcriptomic analysis of the triphasic life history of the model red alga Bostrychia moritziana, as published here: paper ref

## Repository Structure

- `scripts/`: analysis scripts and pipeline runner
- `scripts_by_figure/`: figure- and supplemental-table-level wrappers. Each directory contains a `README.md` and `run.sh`.
- `R/`: shared R helper functions used in the analysis scripts
- `config/`: configuration files and other assets
- `manuscript_outputs/`: final manuscript figures and supplemental tables
- `metadata/`: manifests mapping manuscript outputs to source paths and relevant scripts

## Running Figure Scripts

From the repository root, run a figure-specific script, for example:

```bash
scripts_by_figure/Figure_3/run.sh
```

To run the retained manuscript pipeline in order:

```bash
scripts/run_manuscript_pipeline.sh
```

The scripts use these config variables by default:

```bash
PATHS_CONFIG=config/paths.yml
ANALYSIS_CONFIG=config/analysis.yml
DESIGN_CONFIG=config/design.tsv
CONTRASTS_CONFIG=config/contrasts.tsv
```

The current `config/paths.yml` uses repository-relative placeholder paths for upstream matrices and annotation files and writes analysis outputs to `runs/triphasic_manuscript/results`. Place inputs at those paths or edit the config before rerunning.

## Upstream Processing and Dependencies

This manuscript analysis uses raw STAR/Salmon outputs from the RNA-seq processing pipeline `nf-core/rnaseq` 3.12.0 and Nextflow 23.10.0.

Runtime dependencies, R/Bioconductor package requirements, ASE shell tool requirements, and the exact upstream input paths are documented in `docs/dependencies_and_provenance.md`.

## Reproducibility Status

The repository contains the retained manuscript code, figure/table wrappers, and final copied manuscript outputs. A full source rerun with `scripts/run_manuscript_pipeline.sh` requires the upstream input bundle at the paths configured in `config/paths.yml` and `config/ase.yml`, plus the R/Bioconductor and ASE command-line dependencies listed in `docs/dependencies_and_provenance.md`.

## Manuscript Output Mapping

See `metadata/manuscript_output_manifest.tsv` for the complete mapping from manuscript panel/table to final output file and generating script.
