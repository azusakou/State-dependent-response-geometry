#!/usr/bin/env python3
"""Construct lineage-grounded branch probabilities from LARRY clone descendants."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


CANONICAL = {
    "Undifferentiated": "undifferentiated",
    "Neutrophil": "neutrophil", "Neu": "neutrophil",
    "Monocyte": "monocyte", "Mo": "monocyte",
    "Baso": "basophil", "Ba": "basophil",
    "Mast": "mast", "Ma": "mast",
    "Meg": "megakaryocyte", "Mk": "megakaryocyte",
    "Erythroid": "erythroid", "Er": "erythroid",
    "Eos": "eosinophil", "Ccr7_DC": "dendritic", "cDC": "dendritic",
    "pDC": "plasmacytoid_dc", "pDC1": "plasmacytoid_dc", "pDC2": "plasmacytoid_dc",
    "Lymphoid": "lymphoid",
}


def entropy(probabilities: np.ndarray) -> float:
    values = probabilities[probabilities > 0]
    return float(-(values * np.log(values)).sum())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--clone-matrix", type=Path, required=True)
    parser.add_argument("--clone-targets", type=Path, required=True)
    parser.add_argument("--baseline-assignments", type=Path, required=True)
    parser.add_argument("--terminal-time", type=float, default=6)
    args = parser.parse_args()

    data = ad.read_h5ad(args.h5ad, backed="r")
    obs = data.obs.copy()
    data.file.close()
    clones = sparse.load_npz(args.clone_matrix).tocsc()
    if clones.shape[0] != len(obs):
        raise ValueError(f"clone rows {clones.shape[0]} != cells {len(obs)}")
    obs["canonical_fate"] = obs.cell_type_annotation.astype(str).map(CANONICAL).fillna("other")
    obs["time_numeric"] = pd.to_numeric(obs.time_point, errors="coerce")
    branches = sorted(obs.loc[obs.time_numeric.eq(args.terminal_time), "canonical_fate"].unique())
    terminal_frequency = obs.loc[obs.time_numeric.eq(args.terminal_time), "canonical_fate"].value_counts(normalize=True)
    rare_branches = set(terminal_frequency[terminal_frequency < 0.02].index)

    target_rows = []
    assignment_rows = []
    has_condition = "cytokine_condition" in obs
    for clone_id in range(clones.shape[1]):
        start, end = clones.indptr[clone_id], clones.indptr[clone_id + 1]
        indices = clones.indices[start:end]
        if not len(indices):
            continue
        clone_obs = obs.iloc[indices]
        baseline = clone_obs.loc[clone_obs.time_numeric < args.terminal_time]
        terminal = clone_obs.loc[clone_obs.time_numeric.eq(args.terminal_time)]
        if baseline.empty or terminal.empty:
            continue
        condition_groups = terminal.groupby("cytokine_condition", observed=True) if has_condition else [("standard", terminal)]
        for condition, descendants in condition_groups:
            counts = descendants.canonical_fate.value_counts().reindex(branches, fill_value=0)
            probabilities = counts.to_numpy(float) / counts.sum()
            baseline_count = len(baseline)
            row = {
                "dataset_id": obs.dataset_id.iloc[0],
                "clone_id": clone_id,
                "future_condition": condition,
                "terminal_time": args.terminal_time,
                "baseline_cells": baseline_count,
                "descendant_cells": len(descendants),
                "clone_expansion": len(descendants) / baseline_count,
                "transition_entropy": entropy(probabilities),
                "p_differentiation": 1 - probabilities[branches.index("undifferentiated")] if "undifferentiated" in branches else 1.0,
                "p_rare_fate": float(sum(probabilities[index] for index, branch in enumerate(branches) if branch in rare_branches)),
            }
            for branch, probability in zip(branches, probabilities):
                row[f"p_branch_{branch}"] = probability
            target_rows.append(row)
            target_id = len(target_rows) - 1
            for cell_id, baseline_time in zip(baseline.index.astype(str), baseline.time_numeric):
                assignment_rows.append({
                    "target_id": target_id,
                    "cell_id": cell_id,
                    "clone_id": clone_id,
                    "baseline_time": baseline_time,
                    "future_condition": condition,
                })

    targets = pd.DataFrame(target_rows)
    assignments = pd.DataFrame(assignment_rows)
    args.clone_targets.parent.mkdir(parents=True, exist_ok=True)
    args.baseline_assignments.parent.mkdir(parents=True, exist_ok=True)
    targets.to_csv(args.clone_targets, sep="\t", index_label="target_id")
    assignments.to_csv(args.baseline_assignments, sep="\t", index=False)
    print(
        f"targets={len(targets)} clones={targets.clone_id.nunique()} "
        f"baseline_assignments={len(assignments)} branches={branches} rare={sorted(rare_branches)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
