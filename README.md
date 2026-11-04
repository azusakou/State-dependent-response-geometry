# State-dependent response geometry across public single-cell perturbation resources

Code and compact data for the task catalogue, perturbation-response models, genetic interactions, lineage associations, multiome agreement and gene-level evidence described in the manuscript.

## Run the analysis

Install R with `data.table`, `cluster`, `glmnet`, `nnet`, `igraph`, `xgboost`, `BiocNeighbors`, `R.utils` and `Matrix`. `BiocNeighbors` is available through Bioconductor. Dependency information is in `reference/`.

From the repository root:

```sh
Rscript run_all.R
```

This command refits the models and computes the statistics from the supplied inputs without network access or a GPU. It writes results to `results/`, refreshes the compressed tables and checks them against numerical references. `Rscript run_all.R --quick` skips the 200-draw lineage bootstrap and descendant-depth sensitivity. Both commands evaluate the endpoint-rate network under the 600- and 6,000-iteration settings.

## Repository upload

Extract the archive and upload its contained files and folders together, preserving their relative paths. Each distributed file is below 25 MB and the complete repository contains fewer than 100 files, within the [GitHub browser upload limits](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository).

The Norman data are stored as `inputs/norman_summary.rds` and `inputs/norman_gem_1.rds` through `norman_gem_3.rds`. The summary contains condition metadata and aggregate responses; the three arrays contain technical GEM groups 1–3, 4–6 and 7–8. The input reader assembles them in memory and verifies their checksums, dimensions and labels. No manual joining or large reconstructed file is required.

## Data and tables

`inputs/` contains source-derived cell scores, condition/GEM expression profiles, lineage branch counts and state axes, guide-level RNA/ATAC summaries and fold-specific baseline coordinates. The complete 3,591-gene Norman profiles support interaction reorientation and normalization checks. All 32,738 measured genes are retained in the RNA-baseline profiles so feature selection can be repeated within training folds. Reference objects are used for verification rather than model fitting.

Compressed TSVs in `tables/` can be read without running the analysis:

| Table | Contents |
|---|---|
| `task_catalogue` | Conditions, source/context identifiers, quality gates, exclusions and modality masks. |
| `gene_evidence` | Gene-level evidence, criterion states and the non-targeting-control sensitivity. |
| `interactions` | All 131 Norman pairs, measured single responses, projection effects and technical uncertainty. |
| `mixseq_predictions`, `mixseq_validation_splits` | Predictions, applicability decisions and condition-to-fold assignments. |
| `combination_predictions`, `combination_baseline_predictions` | Gene/family-held-out pair predictions and baseline scores. |
| `multiome_support` | Target/guide counts, agreement statistics and exclusions. |
| `lineage_associations`, `lineage_transfer` | Clone-level permutation inference and cross-resource transfer. |
| `norman_growth` | Guide-matched bulk-growth measurements, RNA basin membership and missingness. |
| `baseline_cell_context_membership` | The 7,328 cells in the 57 primary-task RNA-baseline contexts. |
| `mixseq_nested_predictions`, `mixseq_nested_metrics` | Nested predictions, accepted counts, endpoint errors and matched comparisons. |
| `mixseq_nested_eligibility` | The 120 primary conditions, 144 known-intervention lower-support queries and 66 novel-intervention exclusions. |
| `mixseq_nested_training`, `mixseq_nested_policies` | Training/optimizer audits, selected penalties and uncertainty cutoffs. |

`NA` denotes an unavailable value. Criterion states distinguish an unmeasured assay, insufficient support, an evaluated criterion not met and a met rule. The 715-task primary catalogue and the 705-task non-targeting-control sensitivity are retained separately.

## Count-level scVI training

Use Python 3.12 in an isolated environment:

```sh
python3.12 -m venv .venv-scvi
. .venv-scvi/bin/activate
python -m pip install -r upstream/nested_scvi_requirements.txt
python upstream/train_nested_scvi.py all --bundle . --download --parallel 2 --threads 3
Rscript code/17_nested_scvi_models.R . results/nested_scvi/fits results/scvi_validation
Rscript code/17_nested_scvi_models.R . results/nested_scvi/fits results/scvi_validation_nn6000 --nn-maxit=6000
```

The command downloads five baseline archives, verifies their sizes and MD5 checksums, and selects the 21,936 recorded QC-retained DMSO cells from 169 lines and 335 contexts. The compact cell/gene metadata occupy about 0.6 MB. Source matrices and checkpoints are execution outputs under `results/`.

Five outer cell-line folds each have four inner folds and an outer reference fit: 25 executions with 15 distinct reference-cell subsets. Every held-out line is excluded before variable-gene selection and scVI fitting. The 128 non-benchmark lines form a fixed unlabeled reference panel. Training-only CP10k/log1p values select 4,000 batch-aware Cell Ranger variable genes; scVI uses their UMI counts with 30 latent dimensions, two 128-unit layers, dropout 0.1, gene-wise dispersion and negative-binomial likelihood. Fits use seed 20260713, at most 200 epochs, early-stopping patience 20, validation every five epochs and batch size 256. Query inference uses frozen parameters, and model-state hashes are checked before and after mapping.

The R stage selects scaling, endpoint-specific elastic-net penalties and uncertainty cutoffs within the corresponding training folds. The 600- and 6,000-iteration endpoint-network settings share folds, seeds and other hyperparameters. The 144 lower-support queries use frozen primary-panel models and tuning choices; their outcomes never enter training. The 66 new intervention identities are audited through the predefined applicability rule.

The default R workflow uses `inputs/mixseq_nested_scvi.rds` to refit predictors from the measured fold-specific coordinates. This 375 KB input also retains the selected genes, epoch records and fit audits. The count-level route repeats representation training and supports matching-checkpoint recovery. Numerical differences across CPUs and libraries can affect independently retrained values; the tested runtime is recorded in `reference/`.

## Source reconstruction and figures

Public accessions and source terms are in `upstream/sources.tsv`. `upstream/fetch_public.py` retrieves an explicitly selected source. The MultiPerturb extractor can reconstruct its count subset from the deposited object:

```sh
python3 upstream/fetch_public.py multiperturb
tar -xf downloaded_sources/GSE277747_RAW.tar -C downloaded_sources
Rscript upstream/extract_multiome_counts.R downloaded_sources/GSM8528725_MPS.RDS.gz AUTO rebuilt_counts.rds rebuilt_multiome.rds
```

This extraction needs `SingleCellExperiment` and its dependencies. H5AD extraction scripts also need `HDF5Array`, `DelayedMatrixStats`, `hdf5r` and `Matrix`. After the numerical run, `Rscript upstream/plot_panels.R` draws the combination, duration and evidence panels using `ggplot2` and `patchwork`; `ragg` optionally produces previews.

## Interpretation

The calculations describe conditional response and lineage associations. Transcriptional program tails are not direct death, division or stable-fate measurements. Technical GEM groups are partitions of one pooled experiment. Multiome guide agreement and its negative controls are reported together. The evidence table distinguishes missing measurements from an evaluated criterion not met.

The nested prediction evaluation fixes the archived QC mask, gene availability, program definitions and primary panel. DMSO profiles at the matched horizon are observed query inputs. Source-state clusters and outcome-defined target groups are descriptive partitions. Cluster-bootstrap intervals condition on fixed predictions and decisions, without refitting representations or models. Empirical inner risk tolerances are tuning criteria, not future-risk guarantees. Six-hour/24-hour contrasts sample different populations and do not identify stopping or switching effects.

The bulk-growth assay is a separate ten-day measurement from the Norman study. Its RNA combinations were selected from that screen, so the basin association applies to an ascertained panel. Export identifiers such as `synergistic` and `minimum_local_control_energy` refer to a transcriptomic projection class and a regularized inverse-variance diagnostic.

## Sources and reuse

Credit the investigators listed in `upstream/sources.tsv`. Code is covered by `LICENSE-CODE`; data retain their source-specific terms. The CollecTRI prior incorporates resources with their own reuse conditions.

`inputs/norman_growth_fitness.tsv` retains the LGPLv3 terms of the [Perturbseq_GI source repository](https://github.com/thomasmaxwellnorman/Perturbseq_GI/tree/3b25109aeb9c0c2026bd70abd50304a0ad4e5395). Attribution, source files, fixed commit and the complete license texts are in `LICENSE-NORMAN-GROWTH.txt`.

MIX-Seq-derived inputs retain [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) terms. Cite J. M. McFarland and colleagues, the [MIX-seq dataset](https://doi.org/10.6084/m9.figshare.10298696.v3) and its [source article](https://doi.org/10.1038/s41467-020-17440-w). The cell metadata record the archived QC retention and file-qualified identifiers; reconstruction selects those cells without changing their UMI counts. Source URLs and checksums are in `upstream/mixseq_baseline_sources.tsv`.
