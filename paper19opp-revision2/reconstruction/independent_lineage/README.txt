Independent lineage validation from measured pre-second-treatment clone states

Source
Schaff et al., Clonal differences underlie variable responses to sequential and prolonged treatment, Cell Systems (2024), GSE253739; primary paper https://pmc.ncbi.nlm.nih.gov/articles/PMC11003565/ . The original author Seurat object is referenced in primary_source.sha256 and was read without modification. It contains 74,325 cells and 36,601 gene features. Author Lineage and OG_condition fields define clone identity and treatment history. No Barcode and Still multiple labels are excluded; 39,011 cells have an unambiguous assigned lineage. These are clone units from the source experiment, not independent culture replicates.

Frozen question and analysis
locked_analysis_protocol.txt and its SHA256 were written before new expression scoring or prediction. Counts from the author object, rather than the existing PCA or fate analysis, supply all biological features. All 7,848 normalization genes and all twelve frozen axes are measured. Per-cell log1p CP10k scores are averaged by lineage and treatment. Parent states are measured after the first drug and before the second intervention. The outcome is the later relative's mean apoptosis/proliferation/stress program level under the known second treatment. The primary analysis includes 257 tasks from 166 lineages (parent >= 3 cells, descendant >= 2 cells); strict sensitivity includes 85 tasks from 52 lineages (10/5).

Models
Five outer and four inner ancestral-lineage-grouped folds compare a state-augmented ridge model with treatment/count ridge, treatment-group means and calibrated parent persistence. The best simple comparator and penalties are selected inside training folds. Each lineage has total weight one. Scaling uses training data only. Descendant sampled-cell count is an explicitly post-sampling covariate: this is a conditional measured-response benchmark. Neither missing sampled descendants nor transcriptomic scores establish functional death or a prospective treatment policy.

Results
Primary training-SD-standardized joint RMSE: treatment/count ridge 0.966281; state ridge 0.972738. State-minus-inner-selected-simple difference +0.006458, 95% lineage-bootstrap percentile interval[-0.011114,0.032585]. Strict: 0.991452 versus 1.002056; difference +0.010603 [-0.062194,0.074500]. All 249 bootstrap draws repeat model fitting, scaling, inner penalty selection and simple-comparator selection while retaining complete lineages in fixed grouped folds. Intervals concern the observed captured study population and fold design, not between-experiment replication. Added state features do not demonstrate improved prediction in this analysis; the interval does not establish equivalence or biological irrelevance.

Verification
All 4,334 reconstructed lineage-condition support rows exactly match the earlier metadata-only ledger, verified only after reconstruction. Existing outcome results and PC values were unused. An independent code audit found no material grouping, selection, scaling or mathematical-coupling error. An independent calculation reproduces all ten reported held-out scaled RMSE values within 1e-12; all 2,490 bootstrap model records are finite.

Reproduction
Run inspect_author_object.R, score_and_build_units.R, then fit_grouped_prediction.R from the project directory. The first stage temporarily creates a large raw-count cache; the subsequent stages preserve compact scored cells, clone states, analysis units, grouped folds, predictions, model-selection records and bootstrap outputs. The final analysis can be inspected from analysis_units.rds and the output tables without reopening the large author object. Primary raw data remain at their original source location.

The temporary 1.6-GB raw-count cache was removed after successful scoring and verification; the original author object is unchanged and the cache can be regenerated.
