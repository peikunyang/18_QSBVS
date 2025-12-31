#!/usr/bin/env python3
import os
import argparse

N = 8
BLOCK = N * N * N
TOTAL = 3 * BLOCK

def map_index_ccw(i):
    x = i % N
    y = (i // N) % N
    z = i // (N * N)
    x2 = (N - 1) - y
    y2 = x
    z2 = z
    return z2 * (N * N) + y2 * N + x2

def rotate_block(vals):
    out = [0.0] * BLOCK
    for i in range(BLOCK):
        out[map_index_ccw(i)] = vals[i]
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_root", required=True)
    ap.add_argument("--output_root", required=True)
    args = ap.parse_args()

    os.makedirs(args.output_root, exist_ok=True)

    files = [
        fn for fn in sorted(os.listdir(args.input_root))
        if os.path.isfile(os.path.join(args.input_root, fn))
    ]

    for fn in files:
        in_path = os.path.join(args.input_root, fn)
        rows = []

        with open(in_path, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                z, y, x, v = line.split()
                rows.append((float(z), float(y), float(x), float(v)))

        if len(rows) != TOTAL:
            raise RuntimeError(f"{fn}: expect {TOTAL} rows, got {len(rows)}")

        coords = rows
        values = [r[3] for r in rows]

        out0_vals = values[:]
        out1_vals = []

        for k in range(3):
            block = values[k * BLOCK:(k + 1) * BLOCK]
            out1_vals.extend(rotate_block(block))

        with open(os.path.join(args.output_root, fn + "_0"), "w") as f:
            for (z, y, x, _), v in zip(coords, out0_vals):
                f.write(f"{z:4.1f} {y:4.1f} {x:4.1f} {v:15.12f}\n")

        with open(os.path.join(args.output_root, fn + "_1"), "w") as f:
            for (z, y, x, _), v in zip(coords, out1_vals):
                f.write(f"{z:4.1f} {y:4.1f} {x:4.1f} {v:15.12f}\n")

        print(f"[OK] {fn}")

if __name__ == "__main__":
    main()

