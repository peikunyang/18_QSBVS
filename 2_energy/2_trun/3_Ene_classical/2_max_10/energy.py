#!/usr/bin/env python3
import argparse
import os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--rec_dir", required=True)
ap.add_argument("--lig_dir", required=True)
ap.add_argument("--rec_norm", required=True)
ap.add_argument("--lig_norm", required=True)
ap.add_argument("--outfile", required=True)
args = ap.parse_args()

def load_last_column(path):
    return np.loadtxt(path, dtype=np.float64)[:, -1]

def load_norms(path):
    d = {}
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            k, v = line.split()
            d[k] = float(v)
    return d

rec_norms = load_norms(args.rec_norm)
lig_norms = load_norms(args.lig_norm)

rec_order = ["comp", "apo"]
lig_families = ["PPI-M", "PPI-id", "But1", "But2"]

rows = []

for rec in rec_order:
    rec_path = os.path.join(args.rec_dir, rec)
    if not os.path.isfile(rec_path):
        continue

    r = load_last_column(rec_path)
    rnorm = rec_norms[rec]

    for fam in lig_families:
        lig_files = sorted(
            f for f in os.listdir(args.lig_dir)
            if f.startswith(fam + "_")
        )

        lnorm = lig_norms[fam]

        for lf in lig_files:
            lig_path = os.path.join(args.lig_dir, lf)
            if not os.path.isfile(lig_path):
                continue

            l = load_last_column(lig_path)
            val = float(r @ l) * rnorm * lnorm
            rows.append((rec, lf, val))

with open(args.outfile, "w") as f:
    for rec, lig, val in rows:
        f.write(f"{rec:4s} {lig:10s} {val:20.15f}\n")

print(f"saved {len(rows)} rows to {args.outfile}")

