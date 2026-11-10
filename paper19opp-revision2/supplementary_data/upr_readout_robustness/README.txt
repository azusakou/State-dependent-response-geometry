Curated gene-expression readout robustness

All tables are tab-separated UTF-8; .gz files use gzip compression. NA denotes
an unavailable or mathematically unevaluable quantity, not a biological zero.
table_index.tsv lists every table, its dimensions and purpose.

The fixed panels contain 438 RPE1, 30 Adamson K562, 19 untreated melanoma,
19 IFNg-treated melanoma and 17 co-culture target genes. Genes, not sampled
cells, are the resampling units for association intervals. These datasets and
their original UPR outcomes had already been examined: the new curated sets
provide definition robustness, not new independent validation.

The definitions directory retains all sixteen full, Hallmark-disjoint and unique
readouts, canonical-to-assay memberships, hierarchy/overlap, functional-target
annotation counts and all eight families of one hundred matched random readouts.
Functional RNA annotations and intervention-target membership are distinct;
functional-target strata are sparse, particularly in melanoma. Gene-expression
readouts do not measure canonical pathway activation, ERAD flux or ISR activity.

primary80_rank_associations contains all sixteen readouts in five comparisons,
with signed Spearman estimates, 4,999-gene-bootstrap intervals, 19,999-permutation
raw P values and joint BH80. The four twenty-test subfamilies remain explicit.
secondary160_direction_consistency retains raw signs and source-sign-aligned
external effects; direction nulls preserve both marginal sign frequencies.
Chance-adjusted kappa/MCC intervals are descriptive. The 320 separate minimal
covariate checks retain their complete-case counts and correction family.
Paired definition/context differences and the study-balanced signed-effect
summary are descriptive, without a new significance family or interaction claim.

The primary assay uses the frozen 6,780 measured genes and removes the targeted
gene from both the readout numerator, where applicable, and the normalization
denominator, with target-specific matched-control standardization. The random
readout calibration is a separate fixed-library assay in which BOTH real and
random readouts receive identical numerator exclusion and control standardization.
Its upper-tail fractions are finite-control calibrations, not corrected P values.
All inverse/null readout results and unavailable specificity controls are retained.

The canonical113-gene Hallmark exclusion maps to94 assay symbols after approved-
symbol mapping. The original90-member score/rankings remain unchanged; the added
four-symbol94-member table is an explicitly descriptive identifier sensitivity.
It does not enter the eighty-test family or introduce a new source ranking.

Equal-depth ranges sample target cells without replacement while retaining the
original matched-GEM controls and control SDs. Shared-GEM wild ranges quantify
technical sensitivity conditional on fixed control SDs. For external datasets
without genuine GEM metadata, the table explicitly limits that range to source
technical uncertainty. Neither procedure supplies independent culture replication.

The figure-source directory contains the exact values in the final full-readout,
matched-null and parent-definition panels. Definitions, raw P values, null results
and conditional intervals should be considered together: RPE1 correlations are
not exceptional against matched random readouts, and Hallmark-disjoint Reactome
UPR ordering is not supported. Larger count caches, complete numerical draws,
software and public-source provenance remain in the separate research archive.

The ranking directory contains the retrospective source-quintile gradient and
Precision/NDCG analyses, the joint seven-test Holm family, simultaneous adjacent-
quintile intervals, all fixed null-ranking references and normalization/depth
sensitivities. An overall gradient does not imply monotonic upper-quintile means.
Precision counts signed effects above zero; NDCG uses continuous positive RNA
effect as relevance. Neither metric measures functional benefit. The global87-test
BH table is an additional multiplicity sensitivity retaining the original families.

ranking/alternative_definitions retains source-only alternate-parent lists,
overlap, missing matches, all readout-specific null metrics and paired conditional
comparison intervals. The original Hallmark90 lists remain unchanged; alternate
lists do not replace them. Their memberships/ranking outcomes are not invariant
across definitions. ranking/canonical94_projection applies the original90 source
lists and quintiles to the identifier-corrected94 readout. It introduces no94
source ranking, new P-value family or new null calibration. Filenames referring
to current lists denote these original90 source-frozen lists.

Captured-RNA depth diagnostic

sequencing_depth_covariates gives target-level means of log1p recorded UMI counts
and, separately, log1p the exact common-library count; the schema and availability
tables define their columns. All four existing cell caches contained both fields,
with no missing or invalid values and no additional cell/target filter. These
captured-count proxies differ from the number of sampled cells and from sequenced
read depth or saturation. They can also reflect perturbation-induced RNA content.

sequencing_depth170_conditional_associations retains two separate two-covariate
rank-residual models per comparison: source/validation recorded UMI means, or
source/validation common-library means. It covers the sixteen exact-normalization
curated readouts and the original fixed-library ninety-gene Hallmark reference
in five panels, with a separate joint BH170 family. All earlier test families,
outcomes and rankings remain unchanged. The model-input table contains per-target
response/covariate values, not individual-cell data or numerical draws. Conditional
associations do not isolate pure technical confounding or rescue a failed marginal
or specificity test.

ranking/depth_diagnostic gives descriptive continuous matched gaps after one
all-438-gene recorded-UMI regression with source and RPE1 covariates. Its bootstrap
refits that regression while retaining the original Top10/20 and background lists.
No new P value, adjusted Precision/NDCG or candidate reranking was introduced.
