# Independent perturbation-response validation

This analysis evaluates a discovery-only gene priority against two separately generated perturbation studies. Replogle K562 is the discovery source. Adamson K562 (GSE90546/GSM2406681) is the primary independent-study test; Frangieh melanoma (SCP1064) tests a distinct cellular context under unstimulated, interferon-gamma and co-culture conditions. The three melanoma conditions are separate observations within one study.

All sources use one fixed set of 7,848 measured gene symbols and the exact 12-axis memberships in `final_axis_membership.tsv`. Duplicate raw symbols are summed before computing log1p CP10k. The denominator contains only the shared genes. Three transcriptional programs and nine Hallmark axes are retained. These measurements describe RNA programs, rather than cell death, proliferation events or durable fate changes.

## Reproduce the primary analysis

Use R with `data.table` and `digest`. From this directory:

```sh
Rscript analyze_locked_validation.R .
Rscript run_declared_sensitivities.R .
```

The commands use the supplied frozen discovery ranks, model predictions, comparator vectors and scored external cells. They require no network access, GPU or validation-model fitting. File hashes in `harmonized_lock_hashes.tsv` are checked before analysis. Original ranking and harmonization records are preserved separately from the evaluated common-universe rankings.

The primary score is the absolute dominant discovery-program shift, multiplied by nonnegative agreement across technical GEM partitions and divided by one plus movement on the other two programs. Program identity and direction are determined in the discovery source. Discovery metadata eligibility requires at least 100 retained cells and three GEM groups, without selecting genes by effect significance. Adamson requires at least 100 cells per gene; the primary melanoma panel requires at least 30 cells in the unstimulated condition. These criteria produce 30 and 19 genes, respectively, with 17 retained in the melanoma co-culture comparison. The top 10 are primary; top 20 is evaluated only when a smaller unselected comparison group remains.

Signed external program effects are standardized by matched non-targeting controls. The Adamson primary reference is `63(mod)_pBA580`; `Gal4-4(mod)_pBA582` provides a control-reference sensitivity. The combined-source and pilot-only model priorities are evaluated on the same covered genes: 12 in Adamson and four in melanoma. Missing model predictions are not filled. Full-library normalization, equal-cell resampling, restricted guide assignments and cell-support thresholds are reported as sensitivities. Gene bootstraps condition on the fixed scores and observed panels and do not create biological replication.

Simple comparators include dominant-program magnitude, total response magnitude, technical-group agreement, baseline target expression, discovery cell count and generic essentiality. The essentiality vector is the negative mean available Chronos gene effect across DepMap Public 24Q4 cell lines. It is an external comparator, not another validation endpoint. Unique Entrez-ID mappings are retained; ambiguous or unassayed genes are omitted from comparator-specific matched panels. The dataset is [DepMap 24Q4 Public](https://doi.org/10.25452/figshare.plus.27993248.v1), released under CC BY 4.0; the method is described by [Dempster et al., Genome Biology 2021](https://doi.org/10.1186/s13059-021-02540-7). Source investigators and terms must be credited. Large original matrices are referenced rather than duplicated.

## Reproduce the figure

Use R with `ggplot2`, `patchwork`, `ggrepel`, `svglite` and `ragg` in addition to `data.table`:

```sh
Rscript plot_candidate_validation.R . /path/to/Fig_independent_candidate_validation.pdf
```

The PDF and editable SVG show all primary genes, matched magnitude comparisons, model coverage and the complete 48-axis test family. Panel-specific source tables and the figure legend are included. Raster previews are for layout inspection and need not be distributed.

## Exploratory axis-conservation follow-up

The aggregate-priority test did not demonstrate an external ranking advantage. The subsequent follow-up evaluates all 12 predefined axes across all four external study-condition comparisons. All 48 two-sided gene-permutation tests are jointly BH-adjusted. Per-gene effects, all tests, vector-prediction comparisons, minimal covariate sensitivities, paired context contrasts, leave-one-gene-out estimates and equal-cell analyses are retained. These results are exploratory; they do not retrospectively redefine the primary priority score. Direct response-vector copying is compared with zero and the generic mean of metadata-qualified discovery genes, with no external fitted weights.

To reconstruct the discovery-side 12-axis effects and these follow-ups, retain the companion `cross_dataset` directory containing the scored Replogle cells and source control-expression summaries and run the corresponding scripts with the study root as their argument. The primary ranking evaluation and figure reconstruction above are self-contained within this directory. Raw re-extraction uses the supplied extraction scripts and publicly archived original matrices; the public source accessions and file hashes identify their provenance.

The UPR follow-up is selected from the 48-axis exploratory analysis. Rank conservation is distinct from matching effect direction, transferring response magnitude or establishing a context interaction. All corresponding checks must accompany interpretation. No functional control-node candidate or beneficial intervention is established by a positive RNA-program association.
