Independent HepG2 UPR evaluation

This is the only new holdout: Nadig et al., GSE264667, HepG2. Jurkat outcomes were
not read. Source rankings, the 90-gene readout and evaluation metrics were frozen
before scoring. locked_plan.txt/lock.json preserve the original specification.
The metadata addendum records all 90measured members,6,523/6,780normalizer coverage,
and GBF1/ALG14 cell counts 17/24. Primary>=30-cell list utility remains unestimable;
the>=10-cell analysis was specified from metadata before expression outcomes and
is a lower-support sensitivity, not successful primary full-list validation.

Run in a COPY of the revision directory; retain historical locks/logs separately.
1. Obtain the official matrix with python3 independent_validation/download_parallel.py.
   Internet access is required. The downloader uses two connections with backoff,
   records block SHA256 and verifies cached bytes before resuming. It does not
   assume archived download markers imply a raw file exists. Public whole matrices
   are not duplicated in the archive. The pinned input hash and URL are supplied.
2. python3 independent_validation/run_when_ready.py runs extraction, all frozen
   metrics and shared-GEM sensitivity after download completion. It verifies the
   pinned original input hash. R uses data.table,Matrix,hdf5r,digest; precise R and
   package versions/seeds are retained in actual logs and sessionInfo files.
3. Rscript independent_validation/numerical_QA.R independently reconstructs raw
   count libraries and twenty-four target/mode effects. Its maximum effect difference
   was 6.94e-16; library-count differences were zero.
4. Rscript independent_validation/create_figure.R exports the vector figure and
   exact source tables. On macOS it uses Quartz PDF,svglite and ragg because the
   installed Cairo device lacks X11 libraries. Failed attempts are retained.
5. Rscript independent_validation/export_and_supplement.R reconstructs numerical
   export tables. The delivered index includes an explicitly annotated primary
   placeholder-P field and omits the redundant raw-only interval copy. All raw
   interval bounds remain in the combined raw/projected table and research files.

analyze_holdout.R preserves the sole primary Holm5 family; four statistics are
unestimable and have P=1 only for correction. It reports all three support cohorts
and both normalization modes. Per-gene and cell contrasts are retained. Shared
GEM ranges are conditional technical sensitivities, not culture-level confidence
intervals. Finite matched-module fractions are descriptive, not biological tests.
Frozen source inputs remain under reconstruction/cross_dataset/phase2_upr_rpe1
and reconstruction/upr_rpe1_validation; their historical provenance is preserved.

All source calculations are fixed; no new prediction model is trained. The code
and code hashes predate held-out expression scoring; mechanical plotting/QA and
reproduction repairs are documented without changing scientific choices. Logs
include the original 503 download interruptions, low-concurrency recovery and
successful calculations. 438 genes/90 members yielded876 gene-mode rows. Primary
rank ordering used417 genes. Literal Top 10/20 were not refilled or shortened.


Publication snapshot: scientific calculations and compact results are included; manuscript-editing scripts and correspondence are excluded. The latest source/PDF files are at the revision-directory root. See ../README.md and ../LARGE_INPUTS.tsv for omitted large source-derived inputs.
