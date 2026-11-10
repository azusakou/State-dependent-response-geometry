#!/usr/bin/env python3
"""Evaluate frozen OOF predictions as retrospective transcriptional rankings.

Usage: python3 decision_utility.py INPUT_PREDICTIONS.tsv.gz OUTPUT_DIRECTORY
Only Python's standard library is required. No models or acceptance policies are
refitted, and no outcome selects the candidate sets or resolves prediction ties.
"""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import sys
from collections import defaultdict


def read_tsv(path):
    with gzip.open(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        raise ValueError("No rows to write")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def quantile(values, q):
    v = sorted(values)
    h = (len(v) - 1) * q
    lo, hi = math.floor(h), math.ceil(h)
    return v[lo] + (v[hi] - v[lo]) * (h - lo)


source = Path(sys.argv[1])
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
rows = read_tsv(source)
base = [r for r in rows if r["evaluation_panel"] == "primary"
        and r["nn_iteration_limit"] == "600"]
elastic = {r["pair_id"]: r for r in base
           if r["method"] == "elastic_net"
           and r["selection"] == "original_heuristic_0.15"}
assert len(elastic) == 120
all_contexts = defaultdict(list)
for row in elastic.values():
    all_contexts[row["baseline_context"]].append(row)
multi_contexts = {c: rr for c, rr in all_contexts.items() if len(rr) >= 2}
assert len(all_contexts) == 57 and len(multi_contexts) == 30
policies = [("original_cutoff_0.15", "original_heuristic_0.15", "NA"),
            ("inner_tolerance_0.10", "nested_inner_calibrated", "0.1"),
            ("inner_tolerance_0.15", "nested_inner_calibrated", "0.15"),
            ("inner_tolerance_0.20", "nested_inner_calibrated", "0.2")]
# Fixed desired directions on existing transcriptomic scales; neither is a
# direct measurement of cell death or division, and they are not combined.
endpoints = [("apoptosis", 1), ("proliferation", -1)]
decisions, summary, verification = [], [], []
for policy_name, selection, risk in policies:
    kernel = {r["pair_id"]: r for r in base
              if r["method"] == "kernel_transition"
              and r["selection"] == selection and r["risk_budget"] == risk}
    assert len(kernel) == 120
    eligible = defaultdict(list)
    for pair_id, row in kernel.items():
        enet = elastic[pair_id]
        for endpoint, _ in endpoints:
            assert float(row["observed_" + endpoint]) == float(enet["observed_" + endpoint])
        if row["accepted"] == "TRUE":
            eligible[row["baseline_context"]].append(row)
    candidate_sets = {c: rr for c, rr in eligible.items() if len(rr) >= 2}
    for endpoint, sign in endpoints:
        local = []
        for context, candidates in sorted(candidate_sets.items()):
            assert len({r["perturbation"] for r in candidates}) == len(candidates)
            # Tie-breaking uses recorded intervention identity, never outcomes.
            k_choice = sorted(candidates, key=lambda r: (-sign * float(r["predicted_" + endpoint]), r["perturbation"]))[0]
            e_choice = sorted(candidates, key=lambda r: (-sign * float(elastic[r["pair_id"]]["predicted_" + endpoint]), r["perturbation"]))[0]
            observed = [sign * float(r["observed_" + endpoint]) for r in candidates]
            best = max(observed)
            k_value = sign * float(k_choice["observed_" + endpoint])
            e_value = sign * float(e_choice["observed_" + endpoint])
            local.append({"policy": policy_name, "endpoint": endpoint,
                          "desired_direction": "maximize" if sign == 1 else "minimize",
                          "baseline_context": context, "cell_line": k_choice["cell_line"],
                          "outer_fold": k_choice["outer_fold"],
                          "n_candidates": len(candidates),
                          "candidate_interventions": ";".join(sorted(r["perturbation"] for r in candidates)),
                          "kernel_choice": k_choice["perturbation"],
                          "elastic_choice": e_choice["perturbation"],
                          "kernel_choice_observed_fraction": float(k_choice["observed_" + endpoint]),
                          "elastic_choice_observed_fraction": float(e_choice["observed_" + endpoint]),
                          "kernel_regret": best - k_value,
                          "elastic_regret": best - e_value,
                          "random_expected_regret": best - statistics.mean(observed),
                          "regret_difference_kernel_minus_elastic": e_value - k_value,
                          "kernel_top_choice_hit": int(abs(best-k_value) < 1e-12),
                          "elastic_top_choice_hit": int(abs(best-e_value) < 1e-12),
                          "random_expected_hit": sum(abs(best-v) < 1e-12 for v in observed)/len(observed),
                          "methods_same_choice": int(k_choice["pair_id"] == e_choice["pair_id"])})
        decisions.extend(local)
        clusters = defaultdict(list)
        for item in local:
            clusters[item["cell_line"]].append(item)
        lines = sorted(clusters)
        rng = random.Random(20261106)
        draws = []
        for _ in range(10000):
            resampled = [item for line in rng.choices(lines, k=len(lines)) for item in clusters[line]]
            draws.append(statistics.mean(item["regret_difference_kernel_minus_elastic"] for item in resampled))
        summary.append({"policy": policy_name, "endpoint": endpoint,
                        "n_primary_contexts": len(all_contexts),
                        "n_primary_multi_drug_contexts": len(multi_contexts),
                        "n_matched_candidate_contexts": len(local),
                        "n_cell_lines": len(lines),
                        "n_accepted_tasks": sum(r["accepted"] == "TRUE" for r in kernel.values()),
                        "n_tasks_in_choice_contexts": sum(r["n_candidates"] for r in local),
                        "kernel_top_hits": sum(r["kernel_top_choice_hit"] for r in local),
                        "elastic_top_hits": sum(r["elastic_top_choice_hit"] for r in local),
                        "random_expected_hits": sum(r["random_expected_hit"] for r in local),
                        "same_choice_contexts": sum(r["methods_same_choice"] for r in local),
                        "mean_kernel_regret": statistics.mean(r["kernel_regret"] for r in local),
                        "mean_elastic_regret": statistics.mean(r["elastic_regret"] for r in local),
                        "mean_random_expected_regret": statistics.mean(r["random_expected_regret"] for r in local),
                        "mean_paired_regret_difference": statistics.mean(r["regret_difference_kernel_minus_elastic"] for r in local),
                        "paired_difference_ci_low": quantile(draws, .025),
                        "paired_difference_ci_high": quantile(draws, .975),
                        "n_bootstraps": len(draws),
                        "uncertainty_scope": "cell-line-cluster percentile bootstrap; fixed OOF predictions, candidate sets, endpoints and policies"})
    accepted = [r for r in kernel.values() if r["accepted"] == "TRUE"]
    k_rmse = math.sqrt(statistics.mean(float(r["mse"]) for r in accepted))
    e_rmse = math.sqrt(statistics.mean(float(elastic[r["pair_id"]]["mse"]) for r in accepted))
    verification.append({"policy": policy_name, "accepted": len(accepted),
                         "kernel_rmse_recomputed": k_rmse,
                         "matched_elastic_rmse_recomputed": e_rmse,
                         "difference_recomputed": k_rmse-e_rmse})
write_tsv(out / "decision_contexts.tsv", decisions)
write_tsv(out / "decision_summary.tsv", summary)
write_tsv(out / "rmse_verification.tsv", verification)
protocol = {"input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "seed": 20261106, "bootstrap_draws": 10000,
            "analysis_status": "exploratory post-review analysis of frozen retrospective OOF predictions",
            "candidate_rule": "same recorded baseline_context; at least two kernel-accepted interventions; identical candidates for elastic net",
            "objectives": {"apoptosis": "maximize transcriptional tail fraction", "proliferation": "minimize transcriptional tail fraction"},
            "regret": "best observed desired transcriptomic fraction among accepted candidates minus observed desired fraction for selected intervention",
            "tie_rule": "alphabetic perturbation identity, independent of outcomes",
            "aggregation": "equal weight per eligible context; line-cluster bootstrap preserves multiple contexts per line",
            "interpretation": "does not measure viability benefit, death, division, pretreatment-only policy utility or full candidate-panel utility; candidate support and original QC/program/task definitions remain fixed"}
(out / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
print(json.dumps(summary, indent=2))
