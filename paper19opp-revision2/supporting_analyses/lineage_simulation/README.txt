Lineage sampling simulation

Question
Can descendant sampling change local covariance effective rank when clone-level intrinsic branch probabilities have no stress dependence?

Design
The baseline clone/condition rows, stress and adjustment covariates, descendant counts and stress-excluded PCA neighbor graph are reconstructed from the preserved clone-level source summaries. An R preparation step verifies exact agreement with archived effective ranks and adjusted linear stress coefficients. The Python implementation independently verifies all archived original ranks before simulation.

For each resource and each latent logit standard deviation (0.5, 1, 2), 200 simulations draw one independent standard-normal branch vector per unique clone, sharing that vector across its recorded future conditions. Intrinsic probabilities are softmax(log(pi_condition)+sigma*z_clone). Each pi_condition is the descendant-count-weighted empirical branch frequency with 0.5 added to every branch count. These offsets calibrate branch composition only; the simulated intrinsic probabilities are not estimates of the actual unobserved clone probabilities. No intrinsic dependence on stress or baseline neighborhood is imposed.

Each identical latent draw is observed at: (i) recorded descendant counts; (ii) ten descendants for every existing clone-condition row; and (iii) recorded descendant counts randomly permuted within future condition. An additional scenario replaces baseline neighborhoods by randomly selected condition-matched donors, preserving original neighborhood size and excluding the query. A separate n>=10 population retains the actual 92/1268 or 284/1243 rows, reconstructs PCA/neighborhoods on that selected population and compares intrinsic, recorded-depth and ten-descendant observations. Thus underlying probability vectors remain identical across matched observation scenarios, while counts, donor support and query selection are changed explicitly.

Analysis
For every scenario, local effective rank is exp(entropy(normalized nonnegative covariance eigenvalues)), with exact-zero covariance assigned rank zero as in the manuscript. An optimized Frisch-Waugh regression computes the stress coefficient using residual stress after the manuscript linear nuisance adjustment. Full and selected-population regressions retain their respective original covariate designs, including original log descendant count even in equal-depth scenarios. This choice isolates the outcome-measurement change from changing adjustment. Simulation medians and 2.5th/97.5th percentile ranges describe repeated computational realizations, not biological confidence intervals. No simulation-tail frequency is a formal empirical p value for the biological data: the generative null is deliberately simplified and condition-specific offsets come from the same observed resource.

Interpretation
These controls separate intrinsic probability covariance from finite-sampling geometry and show whether recorded count/neighborhood structure can distort the rank. They do not establish actual intrinsic branch probabilities, infer lost differentiation trajectories, establish a causal effect of stress, or provide independent biological replication. Equal-depth simulations with more cells than observed are hypothetical sensitivity scenarios, not new measured descendants. Selection of deep clones changes the population and its neighborhood geometry; increased coefficient variability in that smaller population is not proof that sampling alone caused the observed sign reversal.

Execution
Rscript prepare_design.R SOURCE_REPRODUCIBILITY_PACKAGE OUTPUT_DIRECTORY
python3 simulate_sampling.py
Rscript plot_simulation.R

The prepared clone-level design, branch probabilities and neighbor tables are included so that the simulation itself does not require raw h5ad matrices. They derive from the preserved source summaries (lineage_inputs.rds), not new biological samples. Metadata records random seed, software version and executed simulation count. Source files outside the revision directory were read only.
