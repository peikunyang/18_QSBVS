#!/usr/bin/env python3
import os
import argparse
import shutil
import sys

N = 8
BLOCK = N * N * N

def read_map(path):
    xs, ys, zs, os_ = [], [], [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            t = line.split()
            if len(t) >= 4:
                xs.append(float(t[0]))
                ys.append(float(t[1]))
                zs.append(float(t[2]))
                os_.append(float(t[3]))
    return xs, ys, zs, os_

def write_map(path, xs, ys, zs, os_):
    with open(path, "w", encoding="utf-8") as w:
        for z, y, x, o in zip(xs, ys, zs, os_):
            w.write(f"{z:8.3f} {y:8.3f} {x:8.3f} {o:14.6f}\n")

def map_index_ccw(i):
    x = i % N
    y = (i // N) % N
    z = i // (N * N)
    x2 = (N - 1) - y
    y2 = x
    z2 = z
    return z2 * (N * N) + y2 * N + x2

def rotate_block_ccw(vals):
    out = [0.0] * BLOCK
    for i in range(BLOCK):
        out[map_index_ccw(i)] = vals[i]
    return out

parser = argparse.ArgumentParser()
parser.add_argument("--input_root", required=True)
parser.add_argument("--output_root", required=True)
args = parser.parse_args()

in_root = os.path.abspath(args.input_root)
out_root = os.path.abspath(args.output_root)

for lig in sorted(os.listdir(in_root)):
    lig_dir = os.path.join(in_root, lig)
    if not os.path.isdir(lig_dir):
        continue

    out_dir0 = os.path.join(out_root, lig + "_0")
    out_dir1 = os.path.join(out_root, lig + "_1")
    os.makedirs(out_dir0, exist_ok=True)
    os.makedirs(out_dir1, exist_ok=True)

    for fn in sorted(os.listdir(lig_dir)):
        in_path = os.path.join(lig_dir, fn)
        if not os.path.isfile(in_path):
            continue

        xs, ys, zs, os_ = read_map(in_path)
        if len(os_) != BLOCK:
            print(f"[ERROR] {in_path}: expect {BLOCK} rows, got {len(os_)}")
            sys.exit(1)

        shutil.copy2(in_path, os.path.join(out_dir0, fn))

        os_rot = rotate_block_ccw(os_)
        write_map(os.path.join(out_dir1, fn), xs, ys, zs, os_rot)

        print(f"[OK] {lig}/{fn} -> {lig}_0/, {lig}_1/")

