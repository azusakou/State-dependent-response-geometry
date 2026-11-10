# Metadata-only search for complete frozen-list support

The primary threshold remains>=30 target cells, with the same90-gene RNA readout, frozen Top10/20 and metrics. No additional expression outcomes were accessed. X-Atlas guide-call counts are upper bounds before RNA QC, not confirmed RNA-cell counts.

| Context | GBF1 | ALG14 | Other blocking support/coverage |
|---|---:|---:|---|
| Jurkat, GSE264667 |8|107|Only88/90 frozen readout genes measured (KDELR3/VEGFA missing). |
| X-Atlas HCT116 |133|229|POLRMT25; RPL7 27 guide-filtered cells. |
| X-Atlas HEK293T |103|341|RPS13 29 guide-filtered cells. |

GARS maps to the author guide-library symbol GARS1, ENSG00000106105. All20 candidate counts, official source metadata, checksums and actual logs are supplied. None of these contexts supports both complete frozen lists under the unchanged primary rule. No contexts were pooled, candidates replaced or standards reduced; the existing HepG2 lower-support sensitivity remains secondary.

From the revision directory, reproduce the X-Atlas metadata check:

```sh
python3 final_adjustment/download_xatlas_metadata.py
python3 final_adjustment/count_xatlas.py
```

The large author guide CSVs are downloaded rather than committed. Their official MD5 values were verified for the original analysis and are in `XAtlas_source_integrity.json`. The source is https://doi.org/10.25452/figshare.plus.29190726 .

For Jurkat:

```sh
python3 final_adjustment/get_jurkat_metadata.py
Rscript final_adjustment/Jurkat_metadata_gate.R
```

Only the first1MiB and last64MiB are fetched into a sparse metadata carrier; range hashes are preserved. This is not a complete expression matrix. Never score missing count bytes. The metadata search does not perform additional hypothesis tests or independent RNA validation.
