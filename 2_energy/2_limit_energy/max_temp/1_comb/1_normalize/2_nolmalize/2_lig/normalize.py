#!/usr/bin/env python3
import argparse, os, sys
import numpy as np

def compute_last_replacement(target_norm, data):
    sum_sq_except = np.sum(data**2) - data[-1]**2
    remaining = target_norm**2 - sum_sq_except
    if remaining < 0:
        return None
    return np.sqrt(remaining)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--norm", type=float, required=True)
    ap.add_argument("--elevate", type=float, required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    T = float(args.norm)
    elev = float(args.elevate)

    for fname in sorted(os.listdir(args.indir)):
        infile = os.path.join(args.indir, fname)
        if not os.path.isfile(infile):
            continue
        try:
            data = np.loadtxt(infile, dtype=np.float64).reshape(-1)
        except Exception as e:
            print(f"[ERROR] load {fname}: {e}", file=sys.stderr)
            sys.exit(1)

        if data.size < 3:
            print(f"[ERROR] {fname}: data length < 3", file=sys.stderr)
            sys.exit(1)

        data = data.copy()
        data[-3] = elev
        r = compute_last_replacement(T, data)
        if r is None:
            print(f"[ERROR] {fname}: no solution for last element", file=sys.stderr)
            sys.exit(1)
        data[-1] = r
        scaled = data / T

        outpath = os.path.join(args.outdir, fname)
        with open(outpath, "w") as f:
            for v in scaled:
                f.write(f"{v:15.12f}\n")

if __name__ == "__main__":
    main()

