#!/usr/bin/env python3
"""
calculate_thresholds.py (Species-Only Version)
Calculates a single global threshold per species across all data.
"""

import os
import csv
import numpy as np
from collections import defaultdict

VALIDATION_CSV = "/rds/general/user/jd1322/home/BirdNET_Validation_Stratified/to_validate.csv"
OUTPUT_DIR     = "/rds/general/user/jd1322/home/BirdNET_Validation_Stratified"
TARGET_PRECISION = 0.9
CONF_STEPS = np.arange(0.1, 1.0, 0.05)


def main():
    print(f"Loading validated detections from {VALIDATION_CSV}...")

    if not os.path.exists(VALIDATION_CSV):
        print(f"ERROR: File does not exist at {VALIDATION_CSV}")
        return

    with open(VALIDATION_CSV, "r", encoding="utf-8-sig", errors="ignore") as f:
        lines = f.readlines()

    if not lines:
        print("ERROR: The file is completely empty.")
        return

    # Autodetect delimiter
    first_line = lines[0]
    delim = ";" if (";" in first_line and "," not in first_line) else ","

    reader = csv.reader(lines, delimiter=delim)
    try:
        header = [h.strip().lower() for h in next(reader)]
    except StopIteration:
        print("ERROR: File contains no header row.")
        return

    try:
        species_idx = header.index("species")
        conf_idx    = header.index("confidence")
        correct_idx = header.index("correct")
    except ValueError:
        print(f"ERROR: Could not find required columns in header: {header}")
        return

    rows = []
    for row in reader:
        if not row or len(row) <= max(species_idx, conf_idx, correct_idx):
            continue

        val = row[correct_idx].strip().lower()
        if val == "":
            continue

        if val in ['1', 'y', 'yes', 'true']:
            correct = 1
        elif val in ['0', 'n', 'no', 'false']:
            correct = 0
        else:
            continue

        try:
            conf = float(row[conf_idx])
            rows.append({
                "species":    row[species_idx].strip(),
                "confidence": conf,
                "correct":    correct,
            })
        except ValueError:
            continue

    print(f"Loaded {len(rows)} validated detections.")

    if not rows:
        print("ERROR: No validated rows could be processed.")
        return

    # Group strictly by SPECIES
    groups = defaultdict(list)
    for row in rows:
        groups[row["species"]].append(row)

    results = []

    for species, dets in sorted(groups.items()):
        n_total = len(dets)
        if n_total < 5:
            print(f"  Skipping {species} — only {n_total} global validated samples")
            continue

        best_threshold = None
        best_precision = None

        for thresh in CONF_STEPS:
            above = [d for d in dets if d["confidence"] >= thresh]
            if len(above) == 0:
                continue
            tp = sum(d["correct"] for d in above)
            precision = tp / len(above)

            if precision >= TARGET_PRECISION and best_threshold is None:
                best_threshold = round(float(thresh), 2)
                best_precision = round(precision, 3)

        baseline = [d for d in dets if d["confidence"] >= 0.1]
        baseline_prec = sum(d["correct"] for d in baseline) / len(baseline) if baseline else 0

        print(f"\nSpecies: {species} (Total Sample Size n={n_total})")
        print(f"  Baseline precision at 0.1: {baseline_prec:.3f}")
        if best_threshold:
            print(f"  Min global threshold for {TARGET_PRECISION} precision: {best_threshold} (precision={best_precision})")
        else:
            print(f"  WARNING: Could not achieve {TARGET_PRECISION} precision with pooled samples")

        results.append({
            "species":             species,
            "n_validated":         n_total,
            "baseline_precision":  round(baseline_prec, 3),
            "threshold_for_0.9":   best_threshold if best_threshold else "not achieved",
            "precision_at_thresh": best_precision if best_precision else "",
        })

    results_path = os.path.join(OUTPUT_DIR, "threshold_results.csv")
    with open(results_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)

    print(f"\nSuccess! Species-only results saved to: {results_path}")


if __name__ == "__main__":
    main()
