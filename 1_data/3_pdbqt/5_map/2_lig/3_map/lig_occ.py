#!/usr/bin/env python3
# Batch ligand charge/occupancy assignment (trilinear) for all files in a directory.
# Output order: x fastest, then y, then z; printed as "z y x value".

import argparse
import os
import numpy as np

def read_ligand(path):
    """Read a single ligand FF file: each line 'type x y z q Rii eps'."""
    types, xyz, q = [], [], []
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(("#", ";", "!")):
                continue
            t = s.split()
            if len(t) < 7:
                continue
            types.append(t[0])
            xyz.append([float(t[1]), float(t[2]), float(t[3])])
            q.append(float(t[4]))
    return np.array(types, dtype=object), np.array(xyz, float), np.array(q, float)

def make_grid(vmin, vmax, spacing):
    xs = np.arange(vmin, vmax + 1e-9, spacing)
    ys = np.arange(vmin, vmax + 1e-9, spacing)
    zs = np.arange(vmin, vmax + 1e-9, spacing)
    return xs, ys, zs

# ---------- WRITE MAP: x fastest, output z y x ----------
def write_map(path, xs, ys, zs, M):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as w:
        for k, z in enumerate(zs):          # z slowest
            for j, y in enumerate(ys):      # y middle
                for i, x in enumerate(xs):  # x fastest
                    w.write(
                        f"{z:8.3f} {y:8.3f} {x:8.3f} {M[i,j,k]:14.6f}\n"
                    )

def bracket_axis(arr, v):
    n = len(arr)
    if n == 1:
        return 0, 0, 0.0
    if v <= arr[0]:
        return 0, 0, 0.0
    if v >= arr[-1]:
        return n-2, n-1, 1.0
    lo, hi = 0, n-1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if arr[mid] <= v:
            lo = mid
        else:
            hi = mid
    denom = (arr[hi] - arr[lo]) if arr[hi] != arr[lo] else 1.0
    t = (v - arr[lo]) / denom
    return lo, hi, float(t)

def add_trilinear(node_grid, xs, ys, zs, x, y, z, val):
    i0,i1,tx = bracket_axis(xs, x)
    j0,j1,ty = bracket_axis(ys, y)
    k0,k1,tz = bracket_axis(zs, z)

    wx0, wx1 = (1.0 - tx), tx
    wy0, wy1 = (1.0 - ty), ty
    wz0, wz1 = (1.0 - tz), tz

    node_grid[i0, j0, k0] += val * (wx0 * wy0 * wz0)
    node_grid[i1, j0, k0] += val * (wx1 * wy0 * wz0)
    node_grid[i0, j1, k0] += val * (wx0 * wy1 * wz0)
    node_grid[i1, j1, k0] += val * (wx1 * wy1 * wz0)
    node_grid[i0, j0, k1] += val * (wx0 * wy0 * wz1)
    node_grid[i1, j0, k1] += val * (wx1 * wy0 * wz1)
    node_grid[i0, j1, k1] += val * (wx0 * wy1 * wz1)
    node_grid[i1, j1, k1] += val * (wx1 * wy1 * wz1)

def derive_out_subdir_name(in_path):
    stem = os.path.splitext(os.path.basename(in_path))[0]
    if stem.endswith("_FF"):
        stem = stem[:-3]
    return stem

def process_one_file(in_path, out_root, xs, ys, zs):
    types, coords, charges = read_ligand(in_path)
    shape = (len(xs), len(ys), len(zs))

    Q    = np.zeros(shape, float)
    N_C  = np.zeros(shape, float)
    N_OA = np.zeros(shape, float)

    if coords.size != 0:
        for (t, (x,y,z), q) in zip(types, coords, charges):
            add_trilinear(Q, xs, ys, zs, x, y, z, q)
            if t == "C":
                add_trilinear(N_C, xs, ys, zs, x, y, z, 1.0)
            elif t == "OA":
                add_trilinear(N_OA, xs, ys, zs, x, y, z, 1.0)

    subdir = derive_out_subdir_name(in_path)
    outdir = os.path.join(out_root, subdir)

    write_map(os.path.join(outdir, "map_charge"), xs, ys, zs, Q)
    write_map(os.path.join(outdir, "map_occ_C"),  xs, ys, zs, N_C)
    write_map(os.path.join(outdir, "map_occ_OA"), xs, ys, zs, N_OA)

    print(f"[OK] {os.path.basename(in_path)} -> {outdir}")

def main():
    ap = argparse.ArgumentParser(
        description="Batch assign ligand charges/occupancies to grid nodes via trilinear weights."
    )
    ap.add_argument("--indir", default=None)
    ap.add_argument("--outroot", default="maps")
    ap.add_argument("--min", type=float, default=-1.2)
    ap.add_argument("--max", type=float, default=1.6)
    ap.add_argument("--spacing", type=float, default=0.4)
    args = ap.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    in_dir = (
        os.path.abspath(os.path.join(script_dir, "..", "..", "1_FF", "FF", "2_lig"))
        if args.indir is None else os.path.abspath(args.indir)
    )

    out_root = os.path.abspath(os.path.join(script_dir, args.outroot))
    os.makedirs(out_root, exist_ok=True)

    if not os.path.isdir(in_dir):
        raise SystemExit(f"[ERROR] Input directory not found: {in_dir}")

    xs, ys, zs = make_grid(args.min, args.max, args.spacing)

    files = [
        os.path.join(in_dir, f)
        for f in sorted(os.listdir(in_dir))
        if not f.startswith(".") and os.path.isfile(os.path.join(in_dir, f))
    ]

    if not files:
        raise SystemExit(f"[WARN] No files found in: {in_dir}")

    for fp in files:
        process_one_file(fp, out_root, xs, ys, zs)

    print("[DONE] All maps generated.")

if __name__ == "__main__":
    main()

