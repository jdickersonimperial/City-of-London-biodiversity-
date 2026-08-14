#!/usr/bin/env python3
"""
validation_stratified_sampler.py

Samples 50 recordings per species across all sites combined,
stratified across confidence bins to cover the full 0.1-1.0 range.
This follows standard validation methodology for BirdNET.

Confidence bins (10 per bin = 50 total):
  0.1-0.3, 0.3-0.5, 0.5-0.7, 0.7-0.9, 0.9-1.0

After manual validation, run calculate_thresholds.py to find
the minimum confidence achieving 0.9 precision per species.
"""

import os
import csv
import glob
import random
import shutil
import argparse
import re
from collections import defaultdict

RESULTS_DIR = "/rds/general/user/jd1322/home/BirdNET_Results_v3"
OUTPUT_DIR  = "/rds/general/user/jd1322/home/BirdNET_Validation_Stratified"
WAV_BASE    = "/rds/general/project/silwood_acoustics/live/CityofLondon"
BASE_2526   = f"{WAV_BASE}/2025-26"

SITE_FOLDERS = {
    "79Charterhouse":  [f"{BASE_2526}/79Charterhouse", f"{WAV_BASE}/Sites02.07.24/79 charterhouse july 2"],
    "120FenchurchSt":  [f"{BASE_2526}/120FenchurchSt", f"{WAV_BASE}/Sites18.06.24/120FenchurchSt"],
    "Aldgate_School":  [f"{BASE_2526}/Aldgate_School"],
    "Barber_Surgeons": [f"{BASE_2526}/Barber_Surgeons", f"{WAV_BASE}/Sites02.07.24/barber surgeons july 2"],
    "Cannon_Bridge":   [f"{BASE_2526}/Cannon_Bridge", f"{WAV_BASE}/Sites09.07.24/Cannon_bridge"],
    "Christchurch":    [f"{BASE_2526}/Christchurch", f"{WAV_BASE}/Sites25.06.24/christchurch data"],
    "Cleary":          [f"{BASE_2526}/Cleary", f"{WAV_BASE}/Sites09.07.24/cleary"],
    "Guildhall":       [f"{BASE_2526}/Guildhall", f"{WAV_BASE}/Sites18.06.24/Guildhall"],
    "Inner_Temple":    [f"{BASE_2526}/Inner_Temple", f"{WAV_BASE}/Sites09.07.24/inner temple"],
    "Mansion_House":   [f"{BASE_2526}/Mansion_House"],
    "Nomura":          [f"{BASE_2526}/Nomura", f"{WAV_BASE}/Sites18.06.24/NomuraGreenRoof"],
    "StDunstan":       [f"{BASE_2526}/StDunstan", f"{WAV_BASE}/Sites18.06.24/StDunstan"],
    "Walbrook_Wharf":  [f"{BASE_2526}/Walbrook_Wharf", f"{WAV_BASE}/Sites18.06.24/WalbrookWharf"],
    "Weil":            [f"{BASE_2526}/Weil", f"{WAV_BASE}/Sites25.06.24/weils data"],
    "Wood_Street":     [f"{BASE_2526}/Wood_Street", f"{WAV_BASE}/Sites09.07.24/wood st"],
}

CONF_BINS  = [(0.1, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]
BIN_LABELS = ["0.1-0.3", "0.3-0.5", "0.5-0.7", "0.7-0.9", "0.9-1.0"]
N_PER_BIN  = 10   # 5 bins x 10 = 50 total per species


def find_wav(filename, site):
    for folder in SITE_FOLDERS.get(site, []):
        matches = glob.glob(os.path.join(folder, "**", filename), recursive=True)
        if matches:
            return matches[0]
    return None


def get_bin(conf):
    for i, (lo, hi) in enumerate(CONF_BINS):
        if lo <= conf < hi:
            return i
    return len(CONF_BINS) - 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-detections", type=int, default=20,
                        help="Min total detections across all sites to include species (default 20)")
    args = parser.parse_args()

    random.seed(42)
    print(f"Stratified validation sampler")
    print(f"  {N_PER_BIN} samples per bin x 5 bins = 50 per species")
    print(f"  Pooled across ALL sites")
    print(f"  Min detections to include species: {args.min_detections}")

    # Load all detections pooled by species (across all sites)
    print("\nLoading detections...")
    # species -> bin_idx -> list of rows
    by_species = defaultdict(lambda: defaultdict(list))

    for csv_file in glob.glob(os.path.join(RESULTS_DIR, "*", "*.csv")):
        try:
            with open(csv_file, newline="") as f:
                for row in csv.DictReader(f):
                    if not row.get("common_name") or not row.get("confidence"):
                        continue
                    try:
                        conf = float(row["confidence"])
                    except:
                        continue
                    if conf < 0.1:
                        continue
                    bin_idx = get_bin(conf)
                    by_species[row["common_name"]][bin_idx].append(row)
        except Exception as e:
            print(f"  WARNING: {csv_file}: {e}")

    print(f"Found {len(by_species)} species")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    validation_rows = []
    total_copied = 0
    missing_wavs = 0

    for species in sorted(by_species.keys()):
        total_dets = sum(len(v) for v in by_species[species].values())
        if total_dets < args.min_detections:
            continue

        species_safe = re.sub(r'[^\w\s-]', '', species).replace(' ', '_')
        print(f"\n{species} (total detections: {total_dets})")

        for bin_idx, label in enumerate(BIN_LABELS):
            bin_dets = by_species[species][bin_idx]
            if not bin_dets:
                print(f"  {label}: no detections")
                continue

            n_sample = min(N_PER_BIN, len(bin_dets))
            sampled  = random.sample(bin_dets, n_sample)

            out_folder = os.path.join(OUTPUT_DIR, species_safe, label)
            os.makedirs(out_folder, exist_ok=True)

            for i, det in enumerate(sampled):
                site     = det["site"]
                wav_path = find_wav(det["file"], site)
                if not wav_path:
                    missing_wavs += 1
                    continue

                conf_str = f"{float(det['confidence']):.3f}"
                out_name = f"{species_safe}_{site}_{det['date']}_conf{conf_str}_{i+1}.wav"
                out_path = os.path.join(out_folder, out_name)

                if not os.path.exists(out_path):
                    shutil.copy2(wav_path, out_path)

                validation_rows.append({
                    "species":    species,
                    "site":       site,
                    "date":       det["date"],
                    "conf_bin":   label,
                    "confidence": det["confidence"],
                    "file":       det["file"],
                    "wav_copy":   out_path,
                    "correct":    "",  # fill in 1=correct, 0=wrong
                })
                total_copied += 1

            print(f"  {label}: {n_sample} samples")

    # Save validation sheet
    val_csv = os.path.join(OUTPUT_DIR, "to_validate.csv")
    if validation_rows:
        with open(val_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=validation_rows[0].keys())
            w.writeheader()
            w.writerows(validation_rows)

    print(f"\nDone!")
    print(f"  Species included: {len(by_species)}")
    print(f"  WAVs copied:      {total_copied}")
    print(f"  Missing WAVs:     {missing_wavs}")
    print(f"  Validation sheet: {val_csv}")
    print(f"\nNext steps:")
    print(f"  1. On your Mac open: /Volumes/RDS/general/user/jd1322/home/BirdNET_Validation_Stratified/")
    print(f"  2. Open to_validate.csv in Excel")
    print(f"  3. Listen to each WAV and fill in 'correct': 1=correct, 0=false positive")
    print(f"  4. Save and run: python ~/scripts/calculate_thresholds.py")


if __name__ == "__main__":
    main()
