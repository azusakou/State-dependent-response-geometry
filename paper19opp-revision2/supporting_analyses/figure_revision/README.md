# Figure reconstruction and specificity analysis

This directory contains the processed inputs and plotting support needed to reproduce the revised source-inventory, benchmark, multiome-specificity, endpoint-graph, duration and gene-evidence figures. It requires no downloads or GPU.

Use R with `data.table`, `ggplot2`, `patchwork`, `Matrix`, `scales` and `svglite`. The tested package versions are recorded in `R_session_info.txt`. From this directory:

```sh
Rscript revise_figures.R
```

The default command reads the bundled `source_inputs` and `plot_support` directories and writes PDF figures, editable SVGs and numerical tables to `reproduced`. To select a source bundle and output directory:

```sh
Rscript revise_figures.R source_inputs /path/to/output_directory
```

The supplied source bundle is the directory containing `inputs`, `reference` and `results`. Plotting support is resolved beside the script, independently of the working directory. PDFs are exported with Quartz on macOS and Cairo elsewhere. Helvetica is requested; equivalent font substitution can change layout on another platform.

The script recomputes the non-targeting-construct comparison directly from archived RNA and ATAC guide responses, applies BH correction across all 27 supported targets, updates task eligibility and gene-evidence classes, and reconstructs the specified figures. Assertions check 10 agreement-selected targets, zero specificity passes, 705 retained tasks and the 389/4,128/5,377 gene counts. The task and gene tables retain historical code identifiers; `causal_regulator` denotes the operational perturbation-supported class and `insufficient_causal_evidence` denotes insufficient perturbation evidence, rather than independently validated causality. Historical classes are explicitly retained in a separate column.

`source_inputs/reference/published_numerical_reference.rds` retains the original processed task, gene and assay-coverage snapshot. It supplies frozen source observations and denominators; the specificity-aware quantities are recomputed by this script. The benchmark, graph and duration TSV inputs are frozen source-level results used to redraw those figures, rather than a rerun of all upstream preprocessing or model training. The complete upstream workflow and original processed-data bundle are available at [State-dependent-response-geometry](https://github.com/azusakou/State-dependent-response-geometry). File checksums identify the exact snapshot used here; no repository commit is inferred.

Source accessions, article DOIs, attribution and available reuse information are retained in `source_inputs/sources.tsv`. MultiPerturb is GSE277747; Norman is GSE133344; MIX-Seq is Figshare 10298696.v3; Replogle is Figshare 20029387; LARRY is GSE140802. The compact RDS files contain derived public-source observations. Their source terms remain applicable. `LICENSE-CODE` covers original analysis code and documentation only and does not relicense third-party data. No Norman bulk-growth raw maps are redistributed in this figure subset.

`SHA256SUMS` covers the reconstruction script and all bundled support/input files. From this directory, verify them with `shasum -a 256 -c SHA256SUMS`. SVGs and numerical tables beside this README are the retained reconstruction outputs; viewing-only raster previews are omitted.
