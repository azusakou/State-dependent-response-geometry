# Source-frozen UPR evaluation in RPE1

This directory contains the statistics for a targeted K562-to-RPE1 cellular-context holdout. RPE1 belongs to the same Replogle study as K562. This follow-up was designed after an exploratory 48-test comparison in separate Adamson and Frangieh studies; it is not an independently preregistered study or a biological-replicate experiment.

The source protocol, all memberships, source ranks, comparator inputs, matches and the primary inference code were frozen before RPE1 expression was read. The 51-file manifest is `../cross_dataset/phase2_upr_rpe1/source_lock_manifest.tsv` (SHA256 `b6d91784f8dd623acfc5e11fdc27730b6f205fed6f8c3448649d10695eb14c52`). The extraction first accessed RPE1 expression at 2026-11-06 04:25:00 UTC, after manifest verification. No primary inference code changed after this lock.

Run from the revision directory, with R and data.table, Matrix and digest installed:

```sh
Rscript reconstruction/upr_rpe1_validation/phase2_statistics_preflight.R
Rscript reconstruction/upr_rpe1_validation/run_phase2_statistics.R reconstruction/cross_dataset/phase2_upr_rpe1 reconstruction/cross_dataset/phase2_upr_rpe1/RPE1_phase2_cache.rds reconstruction/upr_rpe1_validation reconstruction/cross_dataset/phase2_upr_rpe1/source_lock_manifest.tsv
Rscript reconstruction/upr_rpe1_validation/run_phase2_sensitivities.R reconstruction/upr_rpe1_validation
Rscript reconstruction/upr_rpe1_validation/run_equal100_sensitivity.R reconstruction/upr_rpe1_validation
```

The locked runner verifies all source-manifest hashes before loading validation. The synthetic preflight verifies analytic statistics and resampling invariants; its artificial inputs are not biological results. The latter two scripts implement the declared exact-denominator and equal-cell sensitivities and descriptive fixed-list summaries without changing the locked rankings. They were written after primary outcomes became available. Their descriptive sign, absolute-effect and leave-one-set summaries are not new primary hypotheses.

## Measurement and inference

Both contexts use 6,780 measured genes for CP10k and the same 90 UPR members. The own targeted member is removed from the numerator, where applicable, in treated and target-specific control cells. Primary normalization retains the fixed library denominator. The exact-denominator sensitivity recomputes every log term using the library excluding raw target counts, including for targets outside UPR. It does not approximate this change by rescaling an aggregate score.

The primary cohort is the metadata-frozen 438-gene intersection; a stricter RPE1 GEM-support population has 390 genes. Responses are target-cell-weighted matched-GEM contrasts divided by the relevant global control SD. The source score is positive UPR effect multiplied by the fraction of supported source GEMs with a positive contrast. Top10 is nested within the single frozen Top20 list. Five source-covariate-matched backgrounds per selected gene are distinct across all twenty sets.

`primary_Holm3_endpoints.tsv` reports the three primary endpoints with 10,000 two-sided permutations and Holm correction, 2,000 whole-gene/intact-matched-set bootstrap intervals, and conditional shared-GEM wild and crossed intervals. The wild sensitivity uses one Rademacher sign per GEM on weighted, gene-centred observed GEM residuals; the control SD is fixed. Missing gene/GEM observations contribute no residual term. All comparator methods share the same wild weights.

One hundred random 90-gene modules, matched to source-control expression/detection bins and excluding UPR members, remain available in the full module tables. Their own-member numerator correction is identical to real UPR. The primary population has only 16 estimable matched-null priority lists; all 100 correlation tests and all missing-list flags remain. A separately frozen common complete-covariate sensitivity has 370 genes and all 100 matched-null lists. It does not replace the primary population.

Equal100 draws sample 100 target cells per gene without replacement in **both** contexts, 200 times. Original target-specific matched-GEM control means and global control SDs remain fixed. The reported central ranges quantify depth perturbation conditional on those references; they are not gene-bootstrap confidence intervals, independent culture replicates or equal-control-depth analyses. Rankings and matched sets are never refitted.

## Main numerical checks

The primary signed correlation is 0.377505 (gene-bootstrap interval 0.289460–0.461358); the Top20 and Top10 matched signed gaps are 0.328763 and 0.573095. All three Holm-adjusted P values are at most 0.000400. Exact target-denominator exclusion retains these estimates (0.378447, 0.327874, 0.572408), and all 200 equal100 draws retain positive correlations and matched gaps.

The selected lists enrich **relative signed** UPR expression: Top20 mean is 0.044228 versus background −0.284535, with 14/20 positive effects; Top10 mean is 0.243003 versus −0.330092, with 9/10 positive. Descriptive absolute-effect gaps have intervals spanning zero. The positive-source-UPR-effect comparator selects the identical Top10 and has no demonstrated Top20 disadvantage. Agreement weighting therefore adds no established benefit.

UPR correlation is not exceptional relative to the expression/detection-matched random modules. The matched priority utility exceeds all 100 random-module priorities only within the separately frozen 370-gene complete-covariate sensitivity (finite-tail fraction 1/101); the primary panel's 16 evaluable null lists give 1/17. These observations support bounded program-directed relative selection, without calibrated amplitude transport, uniform activation, unique conservation, functional control or beneficial intervention claims.

The companion source directory preserves raw-count caches, source/protocol locks and validation-extraction provenance. Raw public matrices are referenced rather than duplicated. Source-input licenses and terms remain with those companion provenance records; this statistics directory does not relicense them.
