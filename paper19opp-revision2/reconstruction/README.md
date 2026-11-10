# Reproducing the response-validation analyses

This directory preserves the count reconstruction, frozen analysis choices, statistical code and compact outputs supporting the manuscript. It is separate from the minimal LaTeX submission package. Original public matrices are referenced rather than duplicated; `../supplementary_data/table_index.tsv` identifies the final tables accompanying the paper.

## Public inputs

| Source | Public record | Role |
|---|---|---|
| Replogle K562 and RPE1 | [Figshare 20029387](https://doi.org/10.6084/m9.figshare.20029387) | K562 discovery; RPE1 cellular-context holdout from the same study |
| Norman CRISPRa | [GSE133344](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133344) | Independent K562 experiment in the three-study prediction benchmark |
| MIX-Seq | [Figshare version 3](https://doi.org/10.6084/m9.figshare.10298696.v3) | Chemical-response contexts in the three-study benchmark |
| Adamson CRISPRi | [GSE90546 / GSM2406681](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2406681) | Independent-study response validation in K562 |
| Frangieh CRISPR knockout | [SCP1064](https://singlecell.broadinstitute.org/single_cell/study/SCP1064) | Independent melanoma study, with three conditions evaluated separately |
| Schaff lineage tracing | [GSE253739](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253739) | Measured clone state before a second treatment and subsequent lineage-linked RNA profiles |
| DepMap Public 24Q4 | [Figshare dataset](https://doi.org/10.25452/figshare.plus.27993248.v1) | Generic essentiality comparator |

Source filenames, feature annotations, retained cell identities and available checksums are preserved beside each analysis. `cross_dataset/final_source_inventory.tsv`, `external_resources/independent_study_registry.tsv` and `independent_lineage/primary_source.sha256` provide the input map. The registry records the initial resource assessment; the completed Schaff analysis is documented in `independent_lineage/README.txt`. Source licences and repository terms continue to apply.

## Analysis sequence

1. **Common RNA measurements and study-held-out prediction.** `cross_dataset/final_gene_universe.tsv` fixes 7,848 measured gene symbols across the five initial studies. `final_axis_membership.tsv` fixes three signed short modules and nine Hallmark expression summaries. Duplicate raw symbols are summed before log1p CP10k, using the common-gene sum as denominator. `score_common.R` reconstructs scores; `build_analysis_units.R`, `fit_transfer.R`, `paired_transfer_checks.R` and `response_coupling.R` produce the response benchmark and sensitivities. `verify_no_label_leakage.R` checks held-out-outcome isolation. The [component README](cross_dataset/README.md) gives units, normalizations and execution arguments. Files named `harmonized_score_part_*` are the final score reconstruction; preliminary `score_part_*` files are not primary inputs.

2. **Independent-study candidate and axis tests.** `candidate_validation/freeze_harmonized_rank.R` freezes discovery-only priorities. `analyze_locked_validation.R` and `run_declared_sensitivities.R` evaluate the fixed priorities in Adamson and Frangieh. The unsuccessful aggregate-priority tests remain alongside the subsequent complete 48-test axis exploration in `explore_fixed_axis_conservation.R`. `external_resources/` holds Frangieh extraction, metadata and exact own-target numerator/denominator corrections; the corresponding external comparisons are evaluated by `candidate_validation/evaluate_UPR_numerator_exclusion.R` and `evaluate_UPR_exact_composition.R`. See the [candidate-validation README](candidate_validation/README.md).

3. **Source-frozen UPR follow-up in RPE1.** `cross_dataset/phase2_upr_rpe1/` contains the separate 6,780-gene universe, ninety UPR members, source-fixed 438-gene roster, Top10/Top20 lists, comparator covariates, matched backgrounds and 100 expression/detection-matched random modules. Source preparation runs through `01_prepare_source_controls.R`, `02_score_source.R` and `03_prepare_source_ranks.R`. The historical `source_lock_manifest.tsv` and `source_lock_record.txt` preserve the choices frozen before RPE1 expression access. `RPE1_extract_once.R` created the scored validation cache in one sequential pass after verifying that manifest. `upr_rpe1_validation/run_phase2_statistics.R` implements the locked inference; its [README](upr_rpe1_validation/README.md) gives exact commands, dependencies and declared sensitivities. `real_UPR_exact_denominator_sensitivity.R` recomputes real-UPR effects from cached counts using the library excluding each target. `readout_member_sensitivity.R` supplies the explicitly descriptive, fixed-list leave-one-readout-member analysis. The 370-gene complete-covariate sensitivity remains separate from the primary 438-gene family.

4. **Independent measured-lineage prediction.** `independent_lineage/inspect_author_object.R` extracts author counts and metadata; `score_and_build_units.R` reconstructs shared RNA scores and lineage/history units; `fit_grouped_prediction.R` performs nested ancestral-lineage-grouped prediction and full-refitting lineage bootstraps. These scripts run from the parent research-project directory. The temporary raw-count cache can be regenerated; compact scored cells, analysis units, folds, model decisions and final summaries are preserved. The [lineage README](independent_lineage/README.txt) describes the observation design and negative result. `plot_negative_validation.R` generates the supporting figure from these results.

## Re-execution and outputs

The supplied R session-information files record the executed environments. Core packages are `data.table`, `Matrix`, `digest`, `yaml`, `hdf5r` and `HDF5Array`; lineage import additionally uses `SeuratObject`. Figure scripts use `ggplot2`, `patchwork`, `ragg` and `svglite`, with other packages declared in their imports. Downstream statistical reproduction can start from preserved scored caches without reopening the large public matrices.

Raw-extraction scripts record the original read-only input locations and, in some cases, the parent project layout. Configure those locations for a new machine before raw re-extraction. Preserve the historical lock files; relocated or modified scripts should have a separate reproduction manifest rather than overwrite the original record. The source-lock construction script intentionally refuses to recreate a historical pre-outcome lock after validation has been opened.

Final tabular outputs are collected by `prepare_supplementary_data.R`; figure source tables are in `../supplementary_data/figure_source_data/`, and publication PDFs are in `../figure/`. Large caches, diagnostic previews and fitting records belong to this research archive, not the minimal submission ZIP.

RNA-program scores are expression measurements. Relative signed UPR ordering and matched-list enrichment are distinguished from absolute activation, calibrated effect transfer and functional outcomes. The UPR-associated readout does not identify a canonical ER-UPR mechanism. Technical GEM partitions and randomized cell splits are not independent biological replicates. The negative broad-state and initial candidate-priority results are retained with the narrower positive follow-up.
