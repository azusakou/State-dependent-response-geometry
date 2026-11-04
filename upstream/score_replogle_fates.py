#!/usr/bin/env python3
"""Chunked QC and matched-control fate scoring for dense Replogle H5AD files."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml


AXES = ["apoptosis", "proliferation", "stress"]


def mad(values: np.ndarray) -> float:
    center = np.nanmedian(values)
    return float(1.4826 * np.nanmedian(np.abs(values - center)))


def score(matrix: np.ndarray, library: np.ndarray, positions: list[int]) -> np.ndarray:
    if not positions:
        return np.full(len(matrix), np.nan, dtype=np.float32)
    values = np.log1p(matrix[:, positions] / np.maximum(library[:, None], 1) * 1e4)
    return values.mean(axis=1).astype(np.float32)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--axes", type=Path, default=Path("upstream/program_gene_sets.yaml"))
    parser.add_argument("--scdblfinder", type=Path)
    parser.add_argument("--cell-output", type=Path, required=True)
    parser.add_argument("--target-output", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=4000)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--cell-line", default="K562")
    args = parser.parse_args(); specification = yaml.safe_load(args.axes.read_text())["axes"]
    adata = ad.read_h5ad(args.input, backed="r"); obs = adata.obs.copy()
    genes = adata.var.gene_name.astype(str).to_numpy(); gene_index = {gene: index for index, gene in enumerate(genes)}
    positions = {}
    for axis in AXES:
        positive = [gene_index[gene] for gene in specification[axis].get("positive_genes", []) if gene in gene_index]
        negative = [gene_index[gene] for gene in specification[axis].get("negative_genes", []) if gene in gene_index]
        positions[axis] = (positive, negative)
        print(f"{axis}: positive={len(positive)} negative={len(negative)}", flush=True)
    frames = []
    for start in range(0, adata.n_obs, args.chunk_size):
        stop = min(start + args.chunk_size, adata.n_obs)
        matrix = np.asarray(adata.X[start:stop, :], dtype=np.float32)
        library = matrix.sum(axis=1)
        frame = obs.iloc[start:stop][["gem_group", "gene", "gene_id", "transcript", "UMI_count", "mitopercent"]].copy()
        frame.insert(0, "cell_id", obs.index[start:stop].astype(str))
        frame["n_genes_by_counts"] = np.count_nonzero(matrix, axis=1)
        frame["matrix_total_counts"] = library
        for axis, (positive, negative) in positions.items():
            positive_score = score(matrix, library, positive)
            negative_score = score(matrix, library, negative) if negative else np.zeros(len(matrix), dtype=np.float32)
            frame[f"score_{axis}"] = positive_score - negative_score
        frames.append(frame.reset_index(drop=True)); print(f"scored {stop}/{adata.n_obs}", flush=True)
    adata.file.close(); cells = pd.concat(frames, ignore_index=True)
    cells["passes_qc_auto"] = True; threshold_rows = []
    for group, index in cells.groupby("gem_group", observed=True).groups.items():
        index = np.asarray(list(index)); frame = cells.loc[index]
        log_counts = np.log1p(pd.to_numeric(frame.UMI_count, errors="coerce").to_numpy(float))
        log_genes = np.log1p(frame.n_genes_by_counts.to_numpy(float)); mito = frame.mitopercent.to_numpy(float)
        cut = {
            "min_counts": np.expm1(np.median(log_counts) - 3 * mad(log_counts)),
            "max_counts": np.expm1(np.median(log_counts) + 5 * mad(log_counts)),
            "min_genes": np.expm1(np.median(log_genes) - 3 * mad(log_genes)),
            "max_genes": np.expm1(np.median(log_genes) + 5 * mad(log_genes)),
            "max_mitopercent": min(np.quantile(mito, 0.995), np.median(mito) + 5 * mad(mito)),
        }
        passes = (
            frame.UMI_count.between(cut["min_counts"], cut["max_counts"])
            & frame.n_genes_by_counts.between(cut["min_genes"], cut["max_genes"])
            & frame.mitopercent.le(cut["max_mitopercent"])
        )
        cells.loc[index, "passes_qc_auto"] = passes.to_numpy()
        threshold_rows.append({"gem_group": group, **cut})
    if args.scdblfinder and args.scdblfinder.exists():
        calls = pd.read_parquet(args.scdblfinder).set_index("cell_id")
        cells = cells.join(calls[["scDblFinder.score", "scDblFinder.class"]], on="cell_id")
        cells["passes_qc_final"] = cells.passes_qc_auto & cells["scDblFinder.class"].ne("doublet")
    else:
        cells["passes_qc_final"] = cells.passes_qc_auto
    args.cell_output.parent.mkdir(parents=True, exist_ok=True)
    cells.to_parquet(args.cell_output, index=False)
    pd.DataFrame(threshold_rows).to_csv(args.cell_output.with_suffix(".qc_thresholds.tsv"), sep="\t", index=False)
    usable = cells.loc[cells.passes_qc_final].copy(); controls = usable.loc[usable.gene.eq("non-targeting")]
    rows = []
    for (group, perturbation), treated in usable.loc[usable.gene.ne("non-targeting")].groupby(["gem_group", "gene"], observed=True):
        baseline = controls.loc[controls.gem_group.eq(group)]
        if len(treated) < args.min_cells or len(baseline) < args.min_cells:
            continue
        row = {"dataset_id": args.input.stem, "gem_group": group, "cell_line": args.cell_line, "task": "genetic", "baseline": "non-targeting", "perturbation": perturbation, "baseline_cells": len(baseline), "treated_cells": len(treated)}
        for axis in AXES:
            values = baseline[f"score_{axis}"].dropna().to_numpy(float)
            threshold = np.median(values) + mad(values)
            success = int(treated[f"score_{axis}"].gt(threshold).sum()); trials = int(treated[f"score_{axis}"].notna().sum())
            probability = (success + 1) / (trials + 2)
            row[f"{axis}_threshold"] = threshold; row[f"p_{axis}"] = probability
            row[f"p_{axis}_se"] = np.sqrt(probability * (1 - probability) / (trials + 3))
        row["pair_id"] = f"{args.input.stem}|gem{group}|non-targeting->{perturbation}"
        rows.append(row)
    targets = pd.DataFrame(rows); args.target_output.parent.mkdir(parents=True, exist_ok=True)
    targets.to_csv(args.target_output, sep="\t", index=False)
    print(f"cells={len(cells)} final_qc={cells.passes_qc_final.sum()} targets={len(targets)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
