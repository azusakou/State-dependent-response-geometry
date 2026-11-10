UPR quintile and fixed-ranking robustness analyses

These are retrospective extensions of the already evaluated K562-to-RPE1 comparison. RPE1 and K562 are cellular contexts from the same Replogle study. The analyses do not provide new independent validation or functional efficacy evidence. Source groups, ranks, lists and sign thresholds never use RPE1 effects to select a method.

Reproduction
Restore the accompanying reconstruction archive and retain the supplied historical provenance. Run numerical replays in a COPY of the project because result files and run timestamps are written anew. From that copy:

cd paper19opp-revision2
Rscript upr_robustness/ranking/analyze_rank_robustness.R . > upr_robustness/ranking/analysis_stdout.txt 2> upr_robustness/ranking/analysis_stderr.txt
Rscript upr_robustness/ranking/project_fixed_primary_quintiles.R . > upr_robustness/ranking/projection_stdout.txt 2> upr_robustness/ranking/projection_stderr.txt
Rscript upr_robustness/ranking/create_ranking_figure.R . > upr_robustness/ranking/figure_stdout.txt 2> upr_robustness/ranking/figure_stderr.txt
Rscript upr_robustness/ranking/reproduce_alternative_source_rankings.R . > upr_robustness/ranking/source_reproduction_attempt1_stdout.txt 2> upr_robustness/ranking/source_reproduction_attempt1_stderr.txt
Rscript upr_robustness/ranking/evaluate_alternative_rankings.R . > upr_robustness/ranking/alternative_evaluation_attempt1_stdout.txt 2> upr_robustness/ranking/alternative_evaluation_attempt1_stderr.txt
Rscript upr_robustness/ranking/alternative_matched_gap_intervals.R . > upr_robustness/ranking/alternative_matched_interval_stdout.txt 2> upr_robustness/ranking/alternative_matched_interval_stderr.txt
Rscript upr_robustness/ranking/project_canonical94_current_lists.R . > upr_robustness/ranking/canonical94_attempt1_stdout.txt 2> upr_robustness/ranking/canonical94_attempt1_stderr.txt

The source-reproduction wrapper executes the historically locked source calculation into source_reproduction and compares eight source-ranking/matching tables. It leaves the historical alternative_definitions/source_lock_manifest.tsv and first-read records unchanged. Do not rerun the original historical freeze script in its populated output directory. The original freeze command and executed code hashes remain preserved for audit.

Numerical analyses use R; figure export additionally requires svglite and ragg. The observed environment was R 4.4.2 on Apple Silicon; exact package/BLAS versions are in sessionInfo.txt and the corresponding projection/evaluation session records. Seeds are 2026110607 for the original extension and 2026110608 for the definition sensitivities. Draw counts are 2,000 gene/technical resamples, 10,000 permutation/uniform controls, and the original 200 equal-100-cell draws. Plot jitter has a separate seed and cannot alter results.

Definitions and testing
Source-signed-UPR quintiles have fixed sizes 88/88/88/87/87 among 438 genes, with gene-symbol ties. The supplied tables distinguish fixed-primary-group projections from separately labeled within-cohort quintiles. Precision@K is the fraction with signed RPE1 effect d>0; this is a sign criterion, not a minimum meaningful effect or functional phenotype. NDCG@K uses linear relevance max(d,0), discounts 1/log2(rank+1), and the ideal gain for the same cohort as a descriptive normalizer. K is fixed at 10/20.

The sole added inferential family is Holm7: all438 Q5-minus-Q1 mean, Kendall tau-b quintile trend and ordinal slope; on370, priority-minus-raw-UPR Precision@10/@20 and NDCG@10/@20. joint_Holm7_tests.tsv contains raw and adjusted P values. Two-sided permutations retain the original expression-by-count strata and center statistics on their permutation-null mean. Other named-method, random-module, background, definition and coverage comparisons are descriptive; they introduce no new formal P family.

Shared whole-gene weights preserve list overlap. Gene NDCG-difference intervals condition on observed same-cohort ideal gain; shared-GEM draws separately recompute effects, signs and ideal gain. Adjacent-quintile max-t intervals and ordered-draw frequencies distinguish an overall gradient from strict monotonicity. Technical GEM and random-cell resampling do not provide independent culture replication.

Coverage and extensions
The primary roster has 438 genes; fair named-comparator comparisons use the common370 population. The strict390 rule retains only8/10 and17/20 current priorities, so literal TopK utility is NA there: lists are never silently shortened or refilled. Exact-denominator/equal-cell sensitivities retain the original source lists.

alternative_definitions contains source-locked Reactome/GO full and Hallmark-disjoint rankings, unchanged current-list projections, all simple baselines, paired intervals and missing-match coverage. Exact-L-minus-target and fixed-L6780 modes are separate. Proper parent-specific random modules exist only in the fixed mode; exact-mode real results are not calibrated against fixed-mode nulls. canonical94_projection restores four unambiguous identifiers only as a current-list/readout sensitivity, with no94source reranking or null-specificity claim. Original90 results remain unchanged.

Provenance and outputs
Protocols, pre-execution hashes, input/output hashes, actual stdout/stderr, session records, metric tables and computational draws are retained. Source locks preceded the first corresponding alternative-readout evaluation; original RPE1 outcomes were already known. Source/ranking replay checks passed. Failed launch/device attempts are separated in failed_attempts; they changed no scientific choice. Later logs retain package-build-version warnings without numerical errors.

The vector figure Fig_UPR_quintiles_ranking_robustness.pdf is 183×162mm, with embedded Helvetica TrueType text. Editable SVG, native actual-size/300dpi previews, a plain-text legend and source tables accompany it. Visual inspection found no clipping or overlap. Main points: upper quintile means are not strictly ordered; the fixed lists enrich positive RNA responses; agreement weighting has no demonstrated benefit over raw source effects; alternative definitions do not yield invariant candidate priorities.
