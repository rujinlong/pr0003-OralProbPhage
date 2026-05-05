#!/usr/bin/env python3
"""Combine 4 host-prediction methods into a calibrated ensemble.

For each vOTU:
  - collect all (method, host_genus, score) predictions
  - count votes per genus
  - weight by per-method calibration weight (default 1.0; can be tuned via the
    --calibration TSV)
  - top-1 genus by weighted vote share is the ensemble call
  - ensemble_score = winning_genus_weight / sum(all_weights) ∈ [0, 1]
  - n_methods_agreeing = how many methods voted for the winning genus

Inputs (all 4 normalized to: vOTU_id, host_genus, host_score, method):
  --iphop, --spacepharer, --rafah, --blast

Output: ensemble_predictions.tsv with one row per vOTU.

Also writes:
  method_pairwise_overlap.tsv — Jaccard agreement matrix among methods at
                                vOTU-level
  per_vOTU_evidence.tsv       — full long-format evidence per vOTU
"""
import argparse, csv, json
from collections import defaultdict
from itertools import combinations

def read_method(path, method_name):
    if not path or path == "" or path == "NONE":
        return []
    out = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                score = float(row.get("host_score") or "nan")
            except ValueError:
                score = float("nan")
            host = (row.get("host_genus") or "").strip()
            if not host or host.lower() == "unknown":
                continue
            out.append((row["vOTU_id"], host, score, method_name))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iphop");        ap.add_argument("--spacepharer")
    ap.add_argument("--rafah");        ap.add_argument("--blast")
    ap.add_argument("--calibration", help="optional TSV: method, weight (default 1.0 each)")
    ap.add_argument("--out_prefix", default="ensemble")
    args = ap.parse_args()

    weights = defaultdict(lambda: 1.0)
    if args.calibration:
        with open(args.calibration) as fh:
            for line in fh:
                f = line.strip().split()
                if len(f) >= 2:
                    weights[f[0]] = float(f[1])

    records = []
    records += read_method(args.iphop,       "iphop")
    records += read_method(args.spacepharer, "spacepharer")
    records += read_method(args.rafah,       "rafah")
    records += read_method(args.blast,       "blast_imgvr")

    # group by vOTU
    by_votu = defaultdict(list)
    for vid, host, score, method in records:
        by_votu[vid].append((host, score, method))

    # ensemble
    with open(f"{args.out_prefix}_predictions.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["vOTU_id","top_host_genus","ensemble_score",
                    "n_methods_agreeing","contributing_methods","support_json"])
        for vid, evid in sorted(by_votu.items()):
            vote_w = defaultdict(float)
            vote_methods = defaultdict(set)
            total_w = 0.0
            for host, score, method in evid:
                weight = weights[method]
                vote_w[host]      += weight
                vote_methods[host].add(method)
                total_w           += weight
            if total_w == 0: continue
            top_host = max(vote_w, key=vote_w.get)
            score_norm = vote_w[top_host] / total_w
            w.writerow([
                vid, top_host, f"{score_norm:.4f}",
                len(vote_methods[top_host]),
                ",".join(sorted(vote_methods[top_host])),
                json.dumps({h: {"weight": vote_w[h],
                                 "methods": sorted(vote_methods[h])}
                            for h in vote_w})
            ])

    # per-vOTU evidence (long)
    with open(f"{args.out_prefix}_per_vOTU_evidence.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["vOTU_id","method","host_genus","host_score"])
        for vid, host, score, method in records:
            w.writerow([vid, method, host, score])

    # pairwise overlap (Jaccard at vOTU-host level)
    method_set = defaultdict(set)
    for vid, host, score, method in records:
        method_set[method].add((vid, host))
    with open(f"{args.out_prefix}_method_pairwise_overlap.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["method_a","method_b","jaccard","intersect","union"])
        for a, b in combinations(sorted(method_set), 2):
            inter = method_set[a] & method_set[b]
            union = method_set[a] | method_set[b]
            j = len(inter) / len(union) if union else 0.0
            w.writerow([a, b, f"{j:.4f}", len(inter), len(union)])

    print(f"[ensemble] {len(by_votu)} vOTUs scored; outputs: {args.out_prefix}_predictions.tsv etc.")

if __name__ == "__main__":
    main()
