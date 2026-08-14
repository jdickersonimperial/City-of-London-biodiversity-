#!/usr/bin/env python3
"""
all_sites_v3.py
BirdNET analysis for all City of London sites - both 2024 and 2025-26 data.
Each array task processes a chunk of 10 site/date pairs to reduce model loading.
"""

import sys
import os
import csv
import glob
import re
from datetime import datetime

# -------  Site mapping --------
# Each entry: site_name -> list of folder paths to search (both years)

BASE_2526 = "/rds/general/project/silwood_acoustics/live/CityofLondon/2025-26"
BASE_2024 = "/rds/general/project/silwood_acoustics/live/CityofLondon"

SITE_FOLDERS = {
    "79Charterhouse":  [
        f"{BASE_2526}/79Charterhouse",
        f"{BASE_2024}/Sites02.07.24/79 charterhouse july 2",
    ],
    "120FenchurchSt":  [
        f"{BASE_2526}/120FenchurchSt",
        f"{BASE_2024}/Sites18.06.24/120FenchurchSt",
    ],
    "Barber_Surgeons": [
        f"{BASE_2526}/Barber_Surgeons",
        f"{BASE_2024}/Sites02.07.24/barber surgeons july 2",
    ],
    "Cannon_Bridge":   [
        f"{BASE_2526}/Cannon_Bridge",
        f"{BASE_2024}/Sites09.07.24/Cannon_bridge",
    ],
    "Christchurch":    [
        f"{BASE_2526}/Christchurch",
        f"{BASE_2024}/Sites25.06.24/christchurch data",
    ],
    "Cleary":          [
        f"{BASE_2526}/Cleary",
        f"{BASE_2024}/Sites09.07.24/cleary",
    ],
    "Guildhall":       [
        f"{BASE_2526}/Guildhall",
        f"{BASE_2024}/Sites18.06.24/Guildhall",
    ],
    "Inner_Temple":    [
        f"{BASE_2526}/Inner_Temple",
        f"{BASE_2024}/Sites09.07.24/inner temple",
    ],
    "Nomura":          [
        f"{BASE_2526}/Nomura",
        f"{BASE_2024}/Sites18.06.24/NomuraGreenRoof",
    ],
    "StDunstan":       [
        f"{BASE_2526}/StDunstan",
        f"{BASE_2024}/Sites18.06.24/StDunstan",
    ],
    "Walbrook_Wharf":  [
        f"{BASE_2526}/Walbrook_Wharf",
        f"{BASE_2024}/Sites18.06.24/WalbrookWharf",
    ],
    "Weil":            [
        f"{BASE_2526}/Weil",
        f"{BASE_2024}/Sites25.06.24/weils data",
    ],
    "Wood_Street":     [
        f"{BASE_2526}/Wood_Street",
        f"{BASE_2024}/Sites09.07.24/wood st",
    ],
    # 2025-26 only sites (no 2024 equivalent)
    "Aldgate_School":  [f"{BASE_2526}/Aldgate_School"],
    "Mansion_House":   [f"{BASE_2526}/Mansion_House"],
}

OUTPUT_BASE = "/rds/general/user/jd1322/home/BirdNET_Results_v3"
TASKS_FILE  = os.path.join(os.path.expanduser("~"), "all_tasks_v3.csv")
CHUNK_SIZE  = 10
LAT         = 51.5155
LON         = -0.0922
MIN_CONF    = 0.1


def extract_date(filename):
    m = re.search(r'_(\d{8})_', filename)
    return m.group(1) if m else None


def get_all_wavs(site):
    """Get all WAV files across all folders for a site."""
    wavs = []
    for folder in SITE_FOLDERS[site]:
        if os.path.isdir(folder):
            wavs += glob.glob(os.path.join(folder, "**", "*.wav"), recursive=True)
    return wavs


# -------  --list-dates mode --------

if len(sys.argv) == 2 and sys.argv[1] == "--list-dates":
    all_pairs = []
    for site in sorted(SITE_FOLDERS.keys()):
        wavs = get_all_wavs(site)
        if not wavs:
            print(f"Skipping {site} — no WAV files found")
            continue
        dates = sorted(set(
            d for f in wavs
            if (d := extract_date(os.path.basename(f))) is not None
        ))
        print(f"{site}: {len(dates)} dates")
        for d in dates:
            all_pairs.append({"site": site, "date": d})

    # Group into chunks
    chunks = []
    for i in range(0, len(all_pairs), CHUNK_SIZE):
        chunk_id = len(chunks) + 1
        for pair in all_pairs[i:i+CHUNK_SIZE]:
            chunks.append({"chunk_id": chunk_id, "site": pair["site"], "date": pair["date"]})

    n_tasks = chunks[-1]["chunk_id"] if chunks else 0
    print(f"\nTotal site/date pairs: {len(all_pairs)}")
    print(f"Chunk size: {CHUNK_SIZE}")
    print(f"Total tasks: {n_tasks}")

    with open(TASKS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["chunk_id", "site", "date"])
        writer.writeheader()
        writer.writerows(chunks)

    print(f"Saved to {TASKS_FILE}")
    print(f"---> Set #PBS -J 1-{n_tasks} in your .pbs.sh file")
    sys.exit(0)


# -------  Array task mode --------

task_id = int(sys.argv[1])

with open(TASKS_FILE) as f:
    all_rows = list(csv.DictReader(f))

chunk_rows = [r for r in all_rows if int(r["chunk_id"]) == task_id]
if not chunk_rows:
    print(f"No rows for chunk {task_id}")
    sys.exit(0)

print(f"Task {task_id}: {len(chunk_rows)} site/date pairs")

# Check which pairs still need processing
pairs_todo = []
for row in chunk_rows:
    site, date_str = row["site"], row["date"]
    out_csv = os.path.join(OUTPUT_BASE, site, f"{site}_detections_{date_str}.csv")
    if os.path.exists(out_csv):
        print(f"  Already done: {site} {date_str}")
    else:
        pairs_todo.append((site, date_str))

if not pairs_todo:
    print("All done — skipping")
    sys.exit(0)

# Load model once
print("Loading BirdNET model...")
from birdnetlib import Recording
from birdnetlib.analyzer import Analyzer
analyzer = Analyzer()
print("Model loaded!")

for site, date_str in pairs_todo:
    print(f"\n--- {site} {date_str} ---")

    output_path = os.path.join(OUTPUT_BASE, site)
    os.makedirs(output_path, exist_ok=True)
    out_csv = os.path.join(output_path, f"{site}_detections_{date_str}.csv")

    # Find WAV files for this date across all folders
    all_wavs = get_all_wavs(site)
    day_files = [f for f in all_wavs if extract_date(os.path.basename(f)) == date_str]
    print(f"  {len(day_files)} WAV files")

    year  = int(date_str[:4])
    month = int(date_str[4:6])
    day   = int(date_str[6:8])
    rec_date = datetime(year=year, month=month, day=day)

    detections = []
    for f in sorted(day_files):
        try:
            recording = Recording(
                analyzer, f,
                lat=LAT, lon=LON,
                date=rec_date,
                min_conf=MIN_CONF,
            )
            recording.analyze()
            for det in recording.detections:
                detections.append({
                    "site":            site,
                    "date":            date_str,
                    "file":            os.path.basename(f),
                    "start_time":      det["start_time"],
                    "end_time":        det["end_time"],
                    "scientific_name": det["scientific_name"],
                    "common_name":     det["common_name"],
                    "confidence":      det["confidence"],
                })
        except Exception as e:
            print(f"  WARNING {os.path.basename(f)}: {e}")

    fieldnames = ["site","date","file","start_time","end_time",
                  "scientific_name","common_name","confidence"]
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        if detections:
            writer.writerows(detections)

    print(f"  Saved {len(detections)} detections")

print(f"\nChunk {task_id} complete!")
