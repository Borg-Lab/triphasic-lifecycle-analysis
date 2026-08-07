# Gene Master Table Schema

Primary key: `gene_id`

## Core columns
- `gene_id`
- `chr`, `start`, `end`, `strand`
- stage/group mean expression columns (`mean_tpm_*`)
- contrast-level DE fields (`log2fc_*`, `fdr_*`, `base_mean_*`)
- marker fields (`marker_group`, `tau`, `best_log2fc`, `best_fdr`, `marker_flag`)
- module fields (`module_id`, `module_rank`, `module_label`, `kME`, `hub_flag`); prefer dendrogram/order-aware labels for plotting and downstream tables
- annotation fields (`best_hit`, `go_terms`, `interpro` if present)
- evolutionary fields (`phylostratum`, `age_label`, `orthogroup_id`)
- QC fields (`detect_flag`, optional technical flags)
