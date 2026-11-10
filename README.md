# Cross-context transcriptional response ranking

Analysis code and compact inputs for transcriptional response ordering and candidate perturbation prioritization. This repository contains no manuscript or submission files.

## Run

From the repository root:

```sh
# Frozen UPR holdout metrics and technical-GEM sensitivity
Rscript run_upr.R

# Earlier perturbation-response benchmark
Rscript run_all.R
```

The UPR workflow needs R, `data.table` and `digest`. It reproduces correlations, gene-bootstrap intervals, the primary Holm family, Precision@10/20, NDCG, simple baselines,100 matched random-module priorities and shared-GEM sensitivities. Four compact numerical references are checked automatically. Inputs are unchanged; only file paths were adapted from the archived analysis. Outputs go to `results/`.

The earlier workflow additionally needs `cluster`, `glmnet`, `nnet`, `igraph`, `xgboost`, `BiocNeighbors`, `R.utils` and `Matrix`. `Rscript run_all.R --quick` skips its full lineage bootstrap. `upstream/` supplies public-source reconstruction scripts; `tables/` contains compact benchmark outputs.

## Scope and sources

UPR calculations start from released standardized gene/cell effects and frozen source rankings in `inputs/upr/`; raw count reconstruction and the complete research archive are outside this compact repository. The90-gene readout, candidate lists, metrics, seeds and thresholds are preserved. Primary support requires30 target cells. Complete HepG2 Top10/20 utility remains unavailable because GBF1/ALG14 have17/24 cells; the10-cell sensitivity does not replace primary validation. Shared-GEM ranges are technical sensitivities, not culture-replicate confidence intervals. No unique UPR mechanism or advantage over simple source-effect ranking is claimed.

Source studies: [Replogle](https://doi.org/10.6084/m9.figshare.20029387) and [Nadig, GSE264667](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264667). The latter is a separate study with overlapping investigators, not independent-laboratory replication. Other source accessions and terms are in `upstream/sources.tsv`.

Code: `LICENSE-CODE`. Data retain their source terms. Norman growth inputs retain the attribution and LGPLv3 terms in `LICENSE-NORMAN-GROWTH.txt`. MIX-Seq inputs retain [CC BY4.0](https://creativecommons.org/licenses/by/4.0/); cite McFarland and colleagues and the [source dataset](https://doi.org/10.6084/m9.figshare.10298696.v3). Exact file checksums are in `MANIFEST.sha256`.
