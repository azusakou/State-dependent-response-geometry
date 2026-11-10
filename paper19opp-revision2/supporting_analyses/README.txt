Reproducibility archive

This archive contains the additional sampling simulation, retrospective
intervention-choice analysis, specificity-aware evidence reconstruction and
figure sources. It is separate from the minimal LaTeX submission package.

Source workflow and public processed inputs:
https://github.com/azusakou/State-dependent-response-geometry
Source data retain the terms and attribution recorded by that repository.
The preserved complete source bundle remains the input for the optional Norman
profile-level verification; the large expression/GEM arrays are not duplicated
in this additional archive.

The lineage simulation runs directly from its included prepared tables:
  python3 lineage_simulation/simulate_sampling.py
NumPy is required. Its seeds, 200-draw settings and assumptions are in the
script and README.txt. Rebuilding those prepared tables from the full source
bundle uses prepare_design.R with explicit source and output arguments.

The intervention-choice calculation uses only Python's standard library:
  python3 prediction_utility/decision_utility.py prediction_utility/inputs/mixseq_nested_predictions.tsv.gz NEW_OUTPUT
Its fixed out-of-fold inputs are included. The profile-level Norman check uses:
  Rscript prediction_utility/verify_norman.R SOURCE_BUNDLE NEW_OUTPUT
It requires the original norman_summary and three GEM-array input files.
The ranking verification uses SOURCE_BUNDLE/tables and is independent of model
retraining. See prediction_utility/README.md for definitions and restrictions.

Figure reconstruction instructions and compact source inputs are in
figure_revision/README.txt. Simulation figure reconstruction uses:
  Rscript lineage_simulation/plot_simulation.R lineage_simulation NEW_FIGURE.pdf
The executed plotting environment used R, ggplot2, patchwork and Helvetica.

The records verify selected calculations from preserved processed data.
They do not constitute external whole-workflow replication, a new independent
biological study or prospective intervention validation.
