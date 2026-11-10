# Retrospective transcriptional prioritization

This analysis uses the frozen count-level nested scVI out-of-fold predictions in
the preserved processed-data bundle. It tests whether small endpoint-RMSE
differences improve the choice of a measured intervention within the same
cell-line, experiment and horizon. It does not establish prospective therapeutic
utility or improved cell viability.

Run from the manuscript project root, passing the processed-data bundle and an
output directory explicitly:

```sh
python3 decision_utility.py BUNDLE/tables/mixseq_nested_predictions.tsv.gz OUTPUT
Rscript verify_norman.R BUNDLE OUTPUT
python3 verify_combination_ranking.py BUNDLE OUTPUT
```

The Python scripts require only the standard library. The R script uses base R.
It independently assembles the processed condition/GEM response arrays and
reimplements the additive-direction statistics without importing the earlier
analysis functions. This is separate numerical verification from preserved
processed inputs, not independent biological validation or reconstruction of
upstream raw-count preprocessing.

The ranking analysis fixes the original primary panel, trained models, predictions,
acceptance rules and inner-selected policies. Within each recorded baseline
context, at least two distinct kernel-accepted drugs are required. Elastic net
uses exactly the same candidate drugs. Two objectives are evaluated separately:
maximizing the apoptosis-associated transcriptional tail fraction and minimizing
the proliferation-associated tail fraction. Predicted ties are broken by
alphabetical intervention identity, without using outcomes. Regret is the
observed gap between the best accepted drug and the model-selected drug on the
specified transcriptional scale. All eligible contexts have equal weight.

The protocol was defined for this exploratory revision analysis; it was not
preregistered before the original study. The 10,000-draw bootstrap samples cell
lines, preserving their contexts and paired model decisions. Its intervals are
conditional on fixed out-of-fold predictions, candidates and policies; they do not
include encoder training, predictor fitting, policy selection or new-experiment
uncertainty. An exactly zero empirical contrast from identical decisions is not
evidence that future performance is known without uncertainty.

`protocol.json` records the input digest, seed and complete definitions.
`decision_contexts.tsv` contains each decision and its candidate set.
`decision_summary.tsv` contains the aggregate results and conditional intervals.
`rmse_verification.tsv` reproduces the relevant overall endpoint errors directly
from the released predictions. `norman_recomputed_interactions.tsv` and
`norman_verification.tsv` document the separately reimplemented interaction
statistics. `combination_ranking_verification.tsv` recomputes AP and AUROC directly
from the released pair scores.

The original kernel cutoff leaves 17 ranking contexts (62 candidate tasks; 17
lines): 15 contexts from experiment 3 at 24 hours and two from experiment 1 at 24
hours. The 0.10 inner tolerance leaves 16 contexts (56 candidate tasks; 16 lines):
15 from experiment 3 and one from experiment 1 at 24 hours. The original primary
panel has 57 contexts, of which 30 contain at least two measured drugs. Thus these
comparisons concern a restricted common-support choice problem.
