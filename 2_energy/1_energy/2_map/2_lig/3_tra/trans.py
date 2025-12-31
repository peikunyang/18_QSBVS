#!/usr/bin/env python3
import os
import argparse

N = 8
BLOCK = N * N * N

def shift_x_wrap(i):
    x = i % N
    y = (i // N) % N
    z = i // (N * N)
    x2 = (x + 1) % N
    return z * (N * N) + y * N + x2

def read_map(path):
    xs, ys, zs, vs = [], [], [], []
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) == 4:
                xs.append(float(t[0]))
                ys.append(float(t[1]))
                zs.append(float(t[2]))
                vs.append(float(t[3]))
    return xs, ys, zs, vs

def write_map(path, xs, ys, zs, vs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for z, y, x, v in zip(xs, ys, zs, vs):
            f.write(f"{z:4.1f} {y:4.1f} {x:4.1f} {v:15.12f}\n")

def shift_values(vs):
    out = [0.0] * BLOCK
    for i in range(BLOCK):
        out[shift_x_wrap(i)] = vs[i]
    return out

parser = argparse.ArgumentParser()
parser.add_argument("--input_root", required=True)
parser.add_argument("--output_root", required=True)
args = parser.parse_args()

for name in sorted(os.listdir(args.input_root)):
    base = os.path.join(args.input_root, name)
    if not os.path.isdir(base):
        continue

    for tag in ("_0", "_1"):
        os.makedirs(os.path.join(args.output_root, name + tag), exist_ok=True)

    for fn in os.listdir(base):
        in_path = os.path.join(base, fn)
        if not os.path.isfile(in_path):
            continue

        xs, ys, zs, vs = read_map(in_path)
        if len(vs) != BLOCK:
            continue

        write_map(
            os.path.join(args.output_root, name + "_0", fn),
            xs, ys, zs, vs
        )

        write_map(
            os.path.join(args.output_root, name + "_1", fn),
            xs, ys, zs, shift_values(vs)
        )

    print(f"[OK] {name} -> {name}_0 , {name}_1")

