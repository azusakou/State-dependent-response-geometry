#!/usr/bin/env python3
"""Score pathway, regulon and Hallmark axes in an H5AD without densifying X."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--network-dir", type=Path, required=True)
    parser.add_argument("--organism", choices=["human", "mouse"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--obs-columns", default="")
    parser.add_argument("--query", default="")
    parser.add_argument("--counts", action="store_true", help="library-normalize and log1p raw counts before scoring")
    parser.add_argument("--impute-nonfinite", action="store_true", help="replace rare non-finite entries with the per-gene median")
    args = parser.parse_args()

    source = ad.read_h5ad(args.input, backed="r")
    keep = np.ones(source.n_obs, dtype=bool)
    if args.query:
        selected = source.obs.query(args.query).index
        keep = source.obs_names.isin(selected)
    data = source[keep].to_memory()
    source.file.close()
    if args.counts:
        sc.pp.normalize_total(data, target_sum=1e4)
        sc.pp.log1p(data)
    if args.impute_nonfinite:
        if sparse.issparse(data.X):
            if np.any(~np.isfinite(data.X.data)):
                raise ValueError("sparse non-finite imputation is intentionally unsupported; repair the source matrix")
        else:
            matrix = np.asarray(data.X)
            bad = ~np.isfinite(matrix)
            if bad.any():
                medians = np.nanmedian(np.where(np.isfinite(matrix), matrix, np.nan), axis=0)
                medians[~np.isfinite(medians)] = 0.0
                matrix[bad] = medians[np.where(bad)[1]]
                data.X = matrix
                print(f"imputed {bad.sum()} non-finite values with per-gene medians")
    data.var_names = data.var_names.astype(str)
    data.var_names_make_unique()

    score_frames = []
    for resource in ["progeny", "collectri", "hallmark"]:
        net = pd.read_csv(args.network_dir / f"{resource}_{args.organism}.tsv.gz", sep="\t")
        if "weight" not in net:
            net["weight"] = 1.0
        dc.mt.ulm(data, net[["source", "target", "weight"]], tmin=5, raw=False, empty=False, verbose=True)
        score = data.obsm.pop("score_ulm").copy()
        data.obsm.pop("padj_ulm", None)
        score.columns = [f"{resource}__{column}" for column in score.columns]
        score_frames.append(score)
    result = pd.concat(score_frames, axis=1)
    result.index = data.obs_names
    requested = [column for column in args.obs_columns.split(",") if column]
    metadata = data.obs[[column for column in requested if column in data.obs]].copy()
    result = metadata.join(result)
    result.index.name = "cell_id"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.reset_index().to_parquet(args.output, index=False)
    print(f"wrote {args.output}: {result.shape[0]} observations x {len(score_frames)} resource families")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
