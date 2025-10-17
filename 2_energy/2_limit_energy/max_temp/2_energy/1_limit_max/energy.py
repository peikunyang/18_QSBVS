#!/usr/bin/env python3
import argparse, os, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--rec_dir", required=True)
ap.add_argument("--lig_dir", required=True)
ap.add_argument("--outfile", required=True)
args = ap.parse_args()

def list_files(indir):
    return [f for f in sorted(os.listdir(indir)) if os.path.isfile(os.path.join(indir, f))]

rec_files = list_files(args.rec_dir)
lig_files = list_files(args.lig_dir)

rec_vecs = [np.loadtxt(os.path.join(args.rec_dir, f), dtype=np.float64).reshape(-1) for f in rec_files]
lig_vecs = [np.loadtxt(os.path.join(args.lig_dir, f), dtype=np.float64).reshape(-1) for f in lig_files]

rows = []
for i, r in zip(rec_files, rec_vecs):
    for j, l in zip(lig_files, lig_vecs):
        val = float(r @ l)
        rows.append([i, j, val])

with open(args.outfile, "w") as f:
    for i, j, val in rows:
        f.write(f"{i:30s} {j:30s} {val:20.12f}\n")

print(f"saved {len(rows)} rows to {args.outfile}")

