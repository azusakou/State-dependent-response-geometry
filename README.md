# Cross-context ordering and prioritization of UPR-associated transcriptional responses

Code, compact results and manuscript sources for **reproducible ordering of transcriptional program responses across experimental contexts and candidate perturbation prioritization**. The current manuscript and revision analyses are in [`paper19opp-revision2/`](paper19opp-revision2/). The earlier analysis workflow remains available at the repository root.

## Current manuscript and results

- [Manuscript PDF](paper19opp-revision2/manuscript_revised_red.pdf) and [LaTeX source](paper19opp-revision2/manuscript_revised_red.tex).
- [Supplementary PDF](paper19opp-revision2/supplementary_revised_red.pdf), [source](paper19opp-revision2/supplementary_revised_red.tex) and [indexed numerical tables](paper19opp-revision2/supplementary_data/table_index.tsv).
- [Revision analysis guide](paper19opp-revision2/README.md), including public-data access, frozen ranking inputs, independent HepG2 evaluation, baselines and negative controls.

K562 source effects retained signed UPR-associated response ordering in the 438-gene RPE1 holdout (Spearman rho=0.378) and in 417 sufficiently supported perturbations from a separately generated HepG2 screen (rho=0.403; 95% gene interval 0.313–0.489). The HepG2 study has overlapping investigators and is not independent-laboratory replication. Candidate-selection results are conditional on a specified RNA readout. Simple source-UPR-effect ranking remained competitive; neither a uniquely conserved UPR mechanism nor an advantage from agreement weighting was established.

The complete frozen Top10/20 lists **remain unestimable under the primary HepG2 support rule** because GBF1/ALG14 have17/24 cells. The >=10-cell analysis is a separately specified sensitivity and does not replace the >=30-cell primary threshold. A later metadata-only search did not find a context supporting both unchanged lists: Jurkat had8 GBF1 cells and only88/90 readout genes; X-Atlas HCT116 had25 POLRMT/27 RPL7 guide-filtered cells, and HEK293T had29 RPS13 cells. No additional RNA outcomes were read, lists refilled, contexts pooled or standards relaxed.

The former curated-readout main section and figure are retained in Supplementary NoteS13/FigureS15; the independent HepG2 evaluation is main section3.5/Figure4. Manuscript changes remain cumulatively red relative to the final first-revision sources. The manuscript PDFs are the previously compiled and checked17/57-page delivery; this repository update does not change numerical findings or manuscript text.

## Repository contents and verification

Revision2 preserves the existing relative analysis paths and compact frozen inputs. Large public matrices, source-derived caches and large metadata downloads are not committed. [`LARGE_INPUTS.tsv`](paper19opp-revision2/LARGE_INPUTS.tsv) records their exact sizes and SHA256 identities; adjacent reconstruction scripts and source metadata document acquisition. Some downstream analyses require regenerating these inputs or obtaining the full research archive distributed with the manuscript. The repository is not claimed to reproduce every analysis offline without them.

```sh
python3 tools/verify_revision2.py
```

This verifies the distributed revision files against their checksums. For numerical reproduction and compilation, follow the [revision guide](paper19opp-revision2/README.md). Dataset terms and the original investigators' attribution continue to apply; the earlier workflow's source licences are retained below. The `paper19opp-revision1/GitHub` symlink resolves historical source paths to the unchanged legacy inputs/code at this repository root.

## Earlier analysis workflow

The following instructions describe the earlier task-catalogue and response-geometry workflow. They are retained for reproducibility; the current paper's supported contribution and validation scope are given above.

Code and compact data for the task catalogue, perturbation-response models, genetic interactions, lineage associations, multiome agreement and gene-level evidence described in the manuscript.

### Run the analysis

Install R with `data.table`, `cluster`, `glmnet`, `nnet`, `igraph`, `xgboost`, `BiocNeighbors`, `R.utils` and `Matrix`. `BiocNeighbors` is available through Bioconductor. Dependency information is in `reference/`.

From the repository root:

```sh
Rscript run_all.R
```

This command refits the models and computes the statistics from the supplied inputs without network access or a GPU. It writes results to `results/`, refreshes the compressed tables and checks them against numerical references. `Rscript run_all.R --quick` skips the 200-draw lineage bootstrap and descendant-depth sensitivity. Both commands evaluate the endpoint-rate network under the 600- and 6,000-iteration settings.

### Data and tables

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

### Count-level scVI training

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

### Source reconstruction and figures

Public accessions and source terms are in `upstream/sources.tsv`. `upstream/fetch_public.py` retrieves an explicitly selected source. The MultiPerturb extractor can reconstruct its count subset from the deposited object:

```sh
python3 upstream/fetch_public.py multiperturb
tar -xf downloaded_sources/GSE277747_RAW.tar -C downloaded_sources
Rscript upstream/extract_multiome_counts.R downloaded_sources/GSM8528725_MPS.RDS.gz AUTO rebuilt_counts.rds rebuilt_multiome.rds
```

This extraction needs `SingleCellExperiment` and its dependencies. H5AD extraction scripts also need `HDF5Array`, `DelayedMatrixStats`, `hdf5r` and `Matrix`. After the numerical run, `Rscript upstream/plot_panels.R` draws the combination, duration and evidence panels using `ggplot2` and `patchwork`; `ragg` optionally produces previews.

### Interpretation

The calculations describe conditional response and lineage associations. Transcriptional program tails are not direct death, division or stable-fate measurements. Technical GEM groups are partitions of one pooled experiment. Multiome guide agreement and its negative controls are reported together. The evidence table distinguishes missing measurements from an evaluated criterion not met.

The nested prediction evaluation fixes the archived QC mask, gene availability, program definitions and primary panel. DMSO profiles at the matched horizon are observed query inputs. Source-state clusters and outcome-defined target groups are descriptive partitions. Cluster-bootstrap intervals condition on fixed predictions and decisions, without refitting representations or models. Empirical inner risk tolerances are tuning criteria, not future-risk guarantees. Six-hour/24-hour contrasts sample different populations and do not identify stopping or switching effects.

The bulk-growth assay is a separate ten-day measurement from the Norman study. Its RNA combinations were selected from that screen, so the basin association applies to an ascertained panel. Export identifiers such as `synergistic` and `minimum_local_control_energy` refer to a transcriptomic projection class and a regularized inverse-variance diagnostic.

### Sources and reuse

Credit the investigators listed in `upstream/sources.tsv`. Code is covered by `LICENSE-CODE`; data retain their source-specific terms. The CollecTRI prior incorporates resources with their own reuse conditions.

`inputs/norman_growth_fitness.tsv` retains the LGPLv3 terms of the [Perturbseq_GI source repository](https://github.com/thomasmaxwellnorman/Perturbseq_GI/tree/3b25109aeb9c0c2026bd70abd50304a0ad4e5395). Attribution, source files, fixed commit and the complete license texts are in `LICENSE-NORMAN-GROWTH.txt`.

MIX-Seq-derived inputs retain [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) terms. Cite J. M. McFarland and colleagues, the [MIX-seq dataset](https://doi.org/10.6084/m9.figshare.10298696.v3) and its [source article](https://doi.org/10.1038/s41467-020-17440-w). The cell metadata record the archived QC retention and file-qualified identifiers; reconstruction selects those cells without changing their UMI counts. Source URLs and checksums are in `upstream/mixseq_baseline_sources.tsv`.
