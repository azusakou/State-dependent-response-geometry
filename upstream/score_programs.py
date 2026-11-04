#!/usr/bin/env python3
"""Compute compact per-cell fate program scores without modifying count matrices."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


def gene_mean_logcp10k(data: ad.AnnData, genes: list[str]) -> np.ndarray:
    index = data.var_names.get_indexer(genes)
    index = index[index >= 0]
    if not len(index):
        return np.full(data.n_obs, np.nan, dtype=np.float32)
    matrix = sparse.csr_matrix(data.X[:, index], dtype=np.float32)
    library = np.asarray(data.X.sum(axis=1)).ravel().astype(np.float32)
    scale = np.divide(1e4, library, out=np.zeros_like(library), where=library > 0)
    matrix = sparse.diags(scale) @ matrix
    matrix.data = np.log1p(matrix.data)
    return np.asarray(matrix.mean(axis=1)).ravel().astype(np.float32)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--pattern", default="*.qc.h5ad")
    parser.add_argument("--axes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--use-final-qc", action="store_true")
    args = parser.parse_args()
    specification = yaml.safe_load(args.axes.read_text())
    frames = []
    for path in sorted(p for p in args.input_dir.glob(args.pattern) if not p.name.startswith("._")):
        print(f"scoring {path.name}", flush=True)
        data = ad.read_h5ad(path)
        keep = np.ones(data.n_obs, dtype=bool)
        if args.use_final_qc and "passes_qc_final" in data.obs:
            keep = data.obs["passes_qc_final"].fillna(False).astype(bool).to_numpy()
        selected = data[keep]
        columns = [
            column for column in ["dataset_id", "batch", "experiment", "time_hours", "task", "perturbation", "is_control", "cell_line"]
            if column in selected.obs
        ]
        frame = selected.obs[columns].copy()
        frame.insert(0, "cell_id", selected.obs_names.astype(str))
        for axis, details in specification["axes"].items():
            positive = gene_mean_logcp10k(selected, details.get("positive_genes", []))
            negative_genes = details.get("negative_genes", [])
            negative = gene_mean_logcp10k(selected, negative_genes) if negative_genes else np.zeros(selected.n_obs, dtype=np.float32)
            frame[f"score_{axis}"] = positive - negative
        frames.append(frame.reset_index(drop=True))
    result = pd.concat(frames, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.suffix == ".parquet":
        result.to_parquet(args.output, index=False)
    else:
        result.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {len(result)} cells to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
