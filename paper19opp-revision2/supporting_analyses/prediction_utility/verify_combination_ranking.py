#!/usr/bin/env python3
"""Recompute pair-ranking AP and AUROC from released task predictions.

Usage: python3 verify_combination_ranking.py BUNDLE OUTPUT_DIRECTORY
"""
import csv
import gzip
from pathlib import Path
import sys

bundle, output = map(Path, sys.argv[1:])
rows = []
for name in ["combination_predictions", "combination_baseline_predictions"]:
    with gzip.open(bundle / "tables" / (name+".tsv.gz"), "rt") as stream:
        rows.extend(csv.DictReader(stream, delimiter="\t"))
results = []
for split in ["leave_gene", "leave_gene_family"]:
    for method in ["nonlinear_xgboost", "differential_expression_ranking"]:
        data = [r for r in rows if r["split"] == split and r["method"] == method]
        # Older tables spell the family split leave_family.
        if not data and split == "leave_gene_family":
            data = [r for r in rows if r["split"] == "leave_family" and r["method"] == method]
        if not data:
            continue
        data.sort(key=lambda r: -float(r["predicted"]))
        scores = [float(r["predicted"]) for r in data]
        y = [r["observed_binary"] == "TRUE" for r in data]
        positives = sum(y)
        n = len(y)
        tp = 0
        ap = 0
        i = 0
        while i < n:
            j = i + 1
            while j < n and scores[j] == scores[i]:
                j += 1
            group_positive = sum(y[i:j])
            tp += group_positive
            ap += (tp/j) * group_positive/positives
            i = j
        positive_scores = [s for s, label in zip(scores,y) if label]
        negative_scores = [s for s, label in zip(scores,y) if not label]
        auc = sum(float(p>q)+.5*float(p==q) for p in positive_scores for q in negative_scores)/(positives*(n-positives))
        results.append({"split": split, "method": method, "n": n,
                        "positives": positives, "average_precision": ap, "auroc": auc})
with (output/"combination_ranking_verification.tsv").open("w",newline="") as stream:
    writer = csv.DictWriter(stream,fieldnames=list(results[0]),delimiter="\t")
    writer.writeheader();writer.writerows(results)
for row in results:
    print(row)
