#!/usr/bin/env python3
import argparse, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--rec_norm", type=float, required=True)
ap.add_argument("--lig_norm", type=float, required=True)
ap.add_argument("--elevate", type=float, required=True)
ap.add_argument("--rec_file", required=True)
ap.add_argument("--lig_file", required=True)
ap.add_argument("--rec_order", required=True)
ap.add_argument("--lig_order", required=True)
ap.add_argument("--outfile", required=True)
args = ap.parse_args()

def load_vec(path):
    return np.loadtxt(path, dtype=np.float64).reshape(-1)

rec_all = load_vec(args.rec_file)
lig_all = load_vec(args.lig_file)

if rec_all.size != 4096:
    raise ValueError(f"rec size {rec_all.size} != 4096")
if lig_all.size != 32768:
    raise ValueError(f"lig size {lig_all.size} != 32768")

rec = rec_all.reshape(2, 2048)
lig = lig_all.reshape(16, 2048)

with open(args.rec_order) as f:
    rec_names = [line.strip() for line in f if line.strip()]
with open(args.lig_order) as f:
    lig_names = [line.strip() for line in f if line.strip()]

energies = (rec @ lig.T) * args.rec_norm * args.lig_norm * np.sqrt(32.0) - args.elevate

with open(args.outfile, "w") as f:
    for i, rname in enumerate(rec_names):
        for j, lname in enumerate(lig_names):
            f.write(f"{rname:25s} {lname:35s} {energies[i,j]:32.22f}\n")

print(f"saved {2*16} rows to {args.outfile}")

