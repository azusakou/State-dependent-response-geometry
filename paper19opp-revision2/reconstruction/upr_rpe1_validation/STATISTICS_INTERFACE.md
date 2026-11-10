# Phase-2 inference interface

`phase2_statistics.R` contains the fixed statistical functions. It does not read validation expression or select features, gene cohorts, score weights or matching rules.

- `source`: one row per frozen primary gene, with `gene`, `source_delta_z`, `priority_score`, `stratum`. Optional source-only metadata can accompany these fields.
- `validation`: one row per frozen gene, with `gene`, `validation_delta_z`.
- `ranking`: frozen `gene`, integer `rank`, including all primary genes.
- `matches`: five rows per selected gene, with `selected_gene`, `background_gene`; all 120 genes across the 20 six-gene sets must be distinct. Missing matching covariates are support failures, not zero-valued measurements.
- `source_gem` / `validation_gem`: `gene`, `gem`, positive `target_cells`, `contrast_z`. The contrast is the matched-GEM treated-minus-control mean divided by that gene's applicable global control SD. All represented target GEMs contribute with observed-cell weights. No absent gene–GEM combination is treated as an observed zero.
- Comparator rankings: `method`, `gene`, `rank`, restricted identically to the frozen complete-comparator roster. The raw positive-UPR method is identified explicitly.
- Null tables: `module_id`, `gene`, `source_delta_z` or `validation_delta_z`, with null rankings and matched sets retained separately. All 100 modules are reported even if a utility comparison is unestimable.

Primary tests remain the signed-effect Spearman correlation and the frozen Top-20 / Top-10 matched gaps, corrected together by Holm. Association permutations preserve the frozen source expression/count strata. Utility permutations choose the selected role within each frozen six-gene set; the same full set of role draws is used for the nested Top-10 and Top-20 tests.

Whole-gene bootstraps provide association intervals; complete six-gene matched-set bootstraps provide utility intervals. Separate technical sensitivities use a shared Rademacher sign per technical GEM, applied to the observed-cell-weighted per-gene/GEM residuals. The original matched controls and control SD remain fixed. Crossed intervals combine whole-gene or intact matched-set resampling with the same GEM signs. This avoids changing the cohort when a conventional global GEM resample omits every measured partition for a sparse gene. All methods share the same validation GEM draws. These conditional intervals do not establish biological replication or guaranteed finite-sample coverage.

Source ranks, memberships, strata, matches, comparator availability, null memberships and this inference code must be hashed before RPE1 expression is accessed. Any subsequent computational repair must preserve the locked scientific choices and be identified separately.
