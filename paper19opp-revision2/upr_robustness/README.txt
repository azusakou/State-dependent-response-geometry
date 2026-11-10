UPR response ordering: additional robustness analyses

Scope and evidence status
This archive extends the specified current revision only. It does not replace the
study framework, initial discovery analyses or historical source-frozen RPE1
evaluation. All four genetic data sources had already been examined. Definitions,
metrics, source rankings and statistical families were fixed before the respective
new calculations, but the extension is retrospective robustness analysis, not a
fresh independent validation. Original, failed and negative results are retained.

Start here
analysis_scope.txt and scope_lock.json record the scope before new summaries.
current_manuscript_baseline.json and baseline/ identify the exact starting files.
ranking/README.txt, modules/README.md and specificity/README.md give executable
commands, estimands, dependencies, seeds, units and output descriptions. The
historical reconstruction/ directory is supplied alongside this directory in the
same reproducibility archive; do not substitute a different prior revision.

Numerical reproduction
Work in a copy: executed analyses write result tables, logs and new timestamps.
Historical locks must remain unmodified. For replay without downloading the large
public matrices, the supplied scored caches and reconstructed prior cell inputs
support the ranking, inference, resampling and figure commands in the component
READMEs. Full count-level reproduction additionally requires the original public
matrices identified by accession, feature schema and checksum. Configure only
input locations on another machine. No raw matrix was edited or reannotated.

The dependency order is:
1. modules/: authoritative metadata and memberships, source-only expression and
   detection matching, count extraction and score summaries, numerical QA.
2. ranking/: original frozen-list/quintile analysis and source-only alternative
   definition ranks. In the original run these ranks were locked before release
   of corresponding new RPE1 readouts. Replay cannot recreate a historical lock;
   reproduce_alternative_source_rankings.R compares the supplied locked tables.
3. specificity/: primary/secondary inference, paired specification/context checks,
   matched-readout calibration, technical sensitivity and the readout figure.
4. ranking/: alternative-definition evaluation and matched gaps. The separately
   defined modules/canonical_Hallmark94/ identifier sensitivity feeds the two
   descriptive 94-gene projection analyses; it neither changes the original
   90-gene result nor adds a new hypothesis test or source ranking.
5. global_primary_multiplicity.py pools the seven ranking and eighty readout
   primary tests as an additional BH87 sensitivity. It does not replace their
   prespecified Holm7 and BH80 families.
6. prepare_submission_data.R exports compact complete tables and figure source
   data. It adds to the historical Supplementary Data index. Editorial assembly,
   changes and compilation checks are retained separately from scientific results.

Software and execution evidence
Analysis uses mature R/Python tools. Exact R/package/BLAS versions, actual
stdout/stderr, executed-code hashes, source hashes, extraction accounting and
successful numerical checks are retained in the component directories. Failed
launches and mechanical coding repairs are retained and labeled; they did not
change the scientific specifications. Large bootstrap/permutation draws and cell
caches belong here, not in the minimal submission package. Public whole matrices
are not duplicated. Existing GPU workloads were not involved.

Uncertainty and interpretation
Perturbed genes are the primary units; repeated contexts and overlapping lists
are resampled jointly. Cells, GEM groups and computational null draws are not
independent biological replicates. Primary real readouts exclude the targeted
gene from both numerator and library. Matched real/null calibrations use the same
separately labeled fixed-library assay. No exact-library real score is calibrated
against a fixed-library null score. All comparisons retain coverage and missing
values. Literal Top10/20 lists are not filled with replacement genes.

The evidence supports an overall K562-to-RPE1 response gradient and enrichment of
positive effects for the specified RNA readout. It does not establish strict
monotonicity of all five means, uniquely conserved UPR biology, invariant candidate
rankings, agreement-weighting superiority, biochemical activity or functional
intervention benefit. Alternative definitions and negative controls constrain the
claims rather than being discarded. See editorial/claim_evidence_map.tsv and the
revision change record for sentence-to-result traceability.

Submission and research files
The minimal submission ZIP contains TeX sources with embedded references, official
CAS dependencies, used figures, Supplementary Data and the reviewer reply. This
research archive preserves analysis/provenance separately, including the baseline,
definitions, logs, results, figure code and change records. Root PDF files are the
actually compiled review copies. Verification records document compilation from
a freshly extracted submission ZIP as well as PDF geometry and visual checks.

Sequencing-depth addendum
A schema-only check distinguished equal-cell sampling from captured-RNA depth.
The separately locked specificity/analyze_sequencing_depth.R adds 170 conditional
UMI/common-library tests with their own BH170 family, leaving all original test
families unchanged. ranking/depth_diagnostic/analyze_depth_diagnostic.R fits one
locked two-depth nuisance model on the unchanged 438-gene Hallmark outcome and
projects fixed lists/backgrounds onto residuals, with 2,000 whole-gene refits and
no new P, ranking, Precision criterion or causal interpretation. See the respective
addendum locks, actual logs, depth schemas and exported complete tables.
