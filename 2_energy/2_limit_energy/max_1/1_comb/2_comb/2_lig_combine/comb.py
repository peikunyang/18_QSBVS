#!/usr/bin/env python3
import argparse, os, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--indir", required=True)
ap.add_argument("--outfile", required=True)
args = ap.parse_args()

files = sorted(f for f in os.listdir(args.indir) if os.path.isfile(os.path.join(args.indir, f)))
vecs = [np.loadtxt(os.path.join(args.indir, f), dtype=np.float64).reshape(-1) for f in files]
out = np.concatenate(vecs) / np.sqrt(16.0)
np.savetxt(args.outfile, out, fmt="%18.15f")

orderfile = os.path.join(os.path.dirname(args.outfile), "lig_chg_occ_order")
with open(orderfile, "w") as f:
    for fname in files:
        f.write(fname + "\n")

