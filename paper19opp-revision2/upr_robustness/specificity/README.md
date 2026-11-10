# Curated RNA readout robustness

These analyses use previously examined Replogle, Adamson and Frangieh outcomes. New memberships and statistical specifications were fixed before their corresponding summaries; this does not create new independent validation. All original panels and all failed readout definitions remain.

The authoritative source definitions and scoring provenance are in the companion `modules` directory. Four functional GO annotations are retained full, Hallmark-disjoint and unique to each of the four functional sets. Two parent UPR definitions, Reactome and GO, are retained full and Hallmark-disjoint. The sixteen readouts yield eighty source/context comparisons. Their measured sizes are 4–69 genes. They summarize annotated RNA expression and do not measure ERAD/ISR biochemical activity or canonical UPR activation.

Primary measurements retain the frozen 6,780-gene normalizer and exact removal of raw targeted-gene counts from that denominator. An own targeted readout member is also removed from the score numerator in treated and target-specific controls. The matched random-readout assay is a separately labelled fixed-library sensitivity: BOTH real and null scores retain fixed L6780, remove the own readout member, and recompute applicable control SDs. Exact-primary real effects are never compared with fixed-library null effects.

The complete canonical Hallmark set contains 113 genes. Approved-symbol mapping identifies 94 assay symbols; disjoint/null exclusions use this full mapping. The original Hallmark reference remains its frozen 90-member assay and reconstructs the prior source/RPE1 results within numerical precision. A separate 94-member identifier sensitivity is documented by the extraction owner and does not enter the eighty-test family.

## Reproduction

From the revision directory, after the authoritative extraction and source-only ranking locks and the score-release QA:

```sh
Rscript upr_robustness/specificity/prepare_fixed_panels.R .
Rscript upr_robustness/specificity/preflight_statistics.R
Rscript upr_robustness/specificity/analyze_signature_robustness.R . primary
Rscript upr_robustness/specificity/analyze_signature_robustness.R . controls
Rscript upr_robustness/specificity/analyze_specification_controls.R . paired
Rscript upr_robustness/specificity/analyze_specification_controls.R . nulls
Rscript upr_robustness/specificity/analyze_technical_support.R . depth
Rscript upr_robustness/specificity/analyze_technical_support.R . wild
Rscript upr_robustness/specificity/plot_readout_robustness.R .
Rscript upr_robustness/specificity/analyze_Hallmark94_descriptive.R .
```

The primary runner verifies the protocol hash and records response input hashes before loading effects. `protocol_lock.json` and the pre-metric code/definition manifests record timing. Synthetic preflight checks analytic properties without biological files. `coding_repairs.jsonl` documents mechanical repairs without changed scientific choices: double-precision MCC arithmetic, alignment of a raw cache to the already-frozen context-specific panels, and disambiguation of the identifier-sensitivity input axis. Actual successful logs and an analytic reproduction of the integer-overflow failure are retained.

## Tables and interpretation

- `primary80_rank_associations.tsv`: every fixed readout/context rho, conditional 4,999-whole-gene-bootstrap interval, 19,999-draw two-sided gene-permutation P, BH80 and the supplementary within-subfamily BH20. Unavailable tests remain in the fixed correction family. Readout subfamilies and all raw P values remain explicit.
- `secondary160_direction_consistency.tsv` and `chance_adjusted_direction_agreement.tsv`: raw signs, concordance, source-sign-aligned external effect, Cohen kappa and MCC. Direction permutations preserve both sign marginals; no 50% chance assumption is used. Kappa/MCC intervals are descriptive; they add no hypothesis family.
- `all320_minimal_covariate_sensitivities.tsv`: four separately fixed two-covariate rank-residual controls, with complete-case counts and joint BH320. Freedman–Lane tests are conditional diagnostics, not causal adjustment.
- `paired_full_disjoint_unique_specification_differences.tsv` and `paired_Frangieh_context_differences.tsv`: paired whole-gene descriptive intervals. No significant-versus-nonsignificant comparison is used to infer a context interaction.
- `study_balanced_signed_effect_summary.tsv`: a descriptive source-sign-aligned effect summary, averaging melanoma contexts within their study first. It is not an independent-study inverse-variance meta-analysis or calibrated common effect.
- `fixed_library_readout_specificity_calibration.tsv`: all eight available matched-size/expression/detection null families in all five comparisons, including null/inverse results. Finite-null upper-tail fractions are descriptive, not corrected P values. The other functional variants lack a separate size-matched null and are flagged in the availability table.
- `Hallmark94_identifier_sensitivity.tsv`: separately specified, descriptive identifier correction of the original90 assay by four metadata-dictated aliases; paired94-minus90 intervals are not new hypothesis tests or a ranking change.
- `equal_depth_signature_summaries.tsv` and `shared_GEM_technical_signature_intervals.tsv`: fixed-control depth perturbations and conditional technical uncertainty. Equal100 retains all 438 RPE1 targets. External datasets without genuine GEM metadata receive source-GEM-only sensitivity, explicitly labelled as not representing their external technical uncertainty. These ranges do not provide independent culture replication or a coverage guarantee.

The six-full-readout figure, parent full/disjoint intervals and all matched-null calibrations are exported together. The R script, SVG, native preview, actual-PDF preview, legend and three exact source tables are retained. The final PDF embeds editable fonts.

## Findings that bound interpretation

Many curated scores have positive source/RPE1 rank associations, but none of the eight calibrated readouts has exceptional RPE1 correlation relative to its own one hundred matched random readouts. The Reactome UPR Hallmark-disjoint definition fails marginal RPE1 ordering (rho −0.0235, BH80 q 0.760), while the GO disjoint definition is positive. Conditional covariate adjustments do not replace a failed marginal endpoint. The small melanoma panels remain sampling sensitive. Thus these results support definition-dependent expression ordering and challenge a generic, pathway-specific UPR conservation claim. They do not establish a new functional mechanism, broadly beneficial intervention, or definition-independent ranking advantage.

Large cell caches, bootstrap/permutation draws and logs belong in the research archive. Submission data should contain the memberships, fixed panels, complete point/interval/availability tables and exact figure data rather than duplicate raw matrices or previews.
