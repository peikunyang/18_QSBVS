#!/usr/bin/env python3
import sys, os

indir = sys.argv[1]
outfile = sys.argv[2]

target = ("comp", "PPI-M_0_0")

files = [
    "E_max_1",
    "E_max_10",
    "E_max_100",
]

shot_cols = {
    "shot3": 3,
    "shot4": 5,
    "shot5": 7,
    "shot6": 9,
    "shot7": 11,
}

results = []

for fname in files:
    path = os.path.join(indir, fname)

    col_values = {k: [] for k in shot_cols}
    target_vals = {}

    with open(path) as f:
        next(f)
        for line in f:
            if not line.strip():
                continue
            t = line.split()
            key = (t[0], t[1])

            for name, idx in shot_cols.items():
                val = float(t[idx])
                col_values[name].append(val)
                if key == target:
                    target_vals[name] = val

    ranks = {}
    for name in shot_cols:
        sorted_vals = sorted(col_values[name])
        ranks[name] = sorted_vals.index(target_vals[name]) + 1

    results.append((fname, ranks))

with open(outfile, "w") as f:
    f.write(
        f"{'file':10s} "
        + " ".join(f"{k:>6s}" for k in shot_cols)
        + "\n"
    )

    for fname, ranks in results:
        f.write(
            f"{fname:10s} "
            + " ".join(f"{ranks[k]:6d}" for k in shot_cols)
            + "\n"
        )

