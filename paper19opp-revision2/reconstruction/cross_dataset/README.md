# Harmonized response transport across three studies

## Feature lock and source reconstruction

`final_gene_universe.tsv` contains 7,848 unique gene symbols measured in Replogle,
Norman, MIX-Seq, Adamson 10X010 and Frangieh. The intersection uses measurement
metadata only and was fixed before fitting these models or inspecting external
validation effects. Its SHA256 is
`f8924241124f458219e12183948aa7916c86d659bcdb80d4d2460030f5dc884a`.

`final_axis_membership.tsv` freezes three signed short modules and nine Hallmark
gene-expression summaries: UPR, E2F, G2/M, p53, TNF/NF-kB, interferon-gamma,
mTORC1, oxidative phosphorylation and hypoxia. These are mean log-expression
scores, not fitted pathway activities. The membership SHA256 is
`5bad744c5841c10d1863a8061d1865f430e8d650ee8b0c0483aa0c70a0b8d903`.
`final_module_gene_coverage.tsv` includes every specified member and its measured
status. Overlapping gene sets are not independent biological replicates.

Raw counts from duplicated gene symbols are summed before transformation. The
same 7,848-gene count sum defines each cell's CP10k denominator in every source.
Scores are means of log1p CP10k values; negative-module means are subtracted
where specified. Original source QC masks remain fixed. Source paths, sizes,
feature counts and duplicated-symbol counts are in `final_source_inventory.tsv`.
Two Replogle symbols and 95 MIX-Seq symbols are duplicated in their source
annotations and are aggregated. No unmeasured gene is imputed as measured zero.

The final discovery score table has 518,327 cells: 276,108 Replogle, 104,507 Norman
and 137,712 MIX-Seq cells after preferred-control selection. Only the chemical
MIX-Seq branch enters the cross-study prediction benchmark. Cached files named
`harmonized_score_part_*` belong to this final reconstruction. Earlier files named
`score_part_*` and `common_gene_universe.tsv` belong only to the preliminary
three-source normalization sensitivity and are not model inputs.

## Units, targets and leakage controls

The primary benchmark requires at least ten treated cells per intervention and
context and sixty matched controls, without response-magnitude or agreement
filtering. It contains 2,151 Replogle gene/GEM units from 358 genes and 48 GEMs;
2,180 Norman condition/GEM units from 280 conditions and eight GEMs; and 527
MIX-Seq drug/context units from 13 drugs, 169 contexts and 109 cell lines.
Replogle and Norman are separate published K562 experiments. Their GEMs are
technical partitions, not distinct biological contexts. The number of modelled
independent published sources is three, not the number of units, axes or splits.

Each control pool is randomly split into disjoint A, B and C subsets. C supplies
context means and variances. Independently split treated A and B cells, together
with their matching control A and B subsets, supply two response estimates.
The context state is contemporaneous matched-control state, not paired individual
pretreatment state. The reliability target is randomized-cell sampling precision
(squared A/B disagreement), not independent-experiment reproducibility.

`primary_cell_partition_membership.tsv.gz` makes the first split auditable.
`independent_control_parts.tsv` records its control summaries. Ten fixed cell
partitions were analysed; partition 1 is primary and the other nine are
sensitivity analyses. No partition was chosen by its performance.

For each held-out study, feature centring/scaling, study-balanced ridge fitting
and penalty selection use the other two studies only. Penalties are selected
from 0.01, 0.1, 1, 10 and 100 by inner leave-study-out error. Models predicting
response B distinguish C-only state, sampling-support covariates, observed pilot
response A and the incremental C-state features. Thus pilot reconstruction is
explicitly separate from predicting an entirely unmeasured intervention.

`no_label_leakage_checks.tsv` verifies that mutating every held-out B value has
exactly zero effect on prediction or penalty selection, that training and query
unit identifiers are disjoint, and that primary cell partitions do not overlap.

## Interpretation of the completed tests

The attempted common response model did not establish incremental transport
benefit from matched-control state. Three-module RMSE for pilot-only versus
pilot plus state was 0.20205 versus 0.22976 in MIX-Seq, 0.08872 versus 0.09536 in
Norman and 0.14630 versus 0.14698 in Replogle. The unshrunk pilot itself was much
better than transferred ridge in MIX-Seq. The analytic sampling-variance baseline
also outperformed learned disagreement models in log-RMSE across all three
held-out studies. These findings do not establish that biological pretreatment
state is irrelevant: no common intervention-specific input representation was
available, and the genetic control-state variation is largely between technical
lanes.

`primary_response_rmse_by_axis_group.tsv` separates the three core modules, nine
Hallmark scores and all twelve axes. `direct_state_value_cluster_comparisons.tsv`
reports the direct incremental-state contrasts under intervention-only,
context-only and crossed intervention/context bootstrap resampling. MIX-Seq
context resampling keeps complete cell lines together. Norman additionally has
deterministic anchor-gene resampling as a shared-gene sensitivity. These
intervals condition on fixed predictions and cell splits; they do not refit the
pipeline or create biological replication.

`cross_program_coupling.tsv` reports three prespecified independent-half
associations per study, adjusted for context, cell counts and common-library
totals, with BH correction over the nine primary tests. Swapped halves,
alternative cell partitions, leave-context-out intervention magnitude and
Norman anchor-gene clustering are explicitly separate sensitivities. Stress
coupling is not directionally universal across studies. The negative
apoptosis/proliferation relationship is not presented as a novel mechanism.

## Frozen secondary integration test

`integrated_priority_gene_predictions.tsv` exports 358 genes covered by the
held-out-Replogle model trained only on Norman and MIX-Seq. Its frozen SHA256 is
`bb50cfd203e70add1876d7062d01698d4e27237a0dced4f10d45e76c4dda7125`.
Predictions are pooled over eligible Replogle GEMs with treated-cell weights,
then standardized by frozen global harmonized Replogle control SDs. Both
`pilot_plus_state` and its `pilot_response` baseline are exported. The fixed
discovery target, direction and technical agreement remain unchanged. Priority
is the sign-aligned nonnegative target prediction times nonnegative discovery
agreement, divided by one plus the other two predicted programs' norm.
External effects are never read by the exporter. Candidate validation must use
identical covered genes for the integrated score and all comparators; uncovered
genes are not extrapolated or filled. This is a secondary pilot-response test,
not prediction for an intervention with no observed discovery response.

## Re-execution

The preserved `common_cell_scores.rds` and final membership files are sufficient
for downstream reproduction with R and `data.table`. From the revision directory
(with `PROJECT_ROOT` pointing to its parent research project):

```sh
Rscript reconstruction/cross_dataset/build_analysis_units.R PROJECT_ROOT OUTPUT
Rscript reconstruction/cross_dataset/fit_transfer.R PROJECT_ROOT OUTPUT
Rscript reconstruction/cross_dataset/paired_transfer_checks.R PROJECT_ROOT OUTPUT
Rscript reconstruction/cross_dataset/response_coupling.R PROJECT_ROOT OUTPUT
Rscript reconstruction/cross_dataset/export_integrated_priority_inputs.R PROJECT_ROOT OUTPUT
Rscript reconstruction/cross_dataset/verify_no_label_leakage.R PROJECT_ROOT OUTPUT
```

The raw reconstruction additionally needs `yaml`, `Matrix`, `HDF5Array` and
`hdf5r`, the source count files identified in the manifest, the external measured
gene metadata and the verified local Hallmark snapshot. `score_common.R` records
its read-only source locations. No raw source matrix is modified. Analysis
outputs are separate from the minimal LaTeX submission package.
