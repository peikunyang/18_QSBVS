#!/usr/bin/env python3
import os
import sys
import math

root = sys.argv[1]

out_dir = "maps"
os.makedirs(out_dir, exist_ok=True)

names = ["map_charge", "map_occ_C", "map_occ_OA"]
lengths = []

for sub in ["PPI-M", "PPI-id", "But1", "But2"]:
    base = os.path.join(root, sub)

    coords = []
    values = []

    with open(os.path.join(base, names[0])) as f:
        for line in f:
            x, y, z, v = line.split()
            coords.append((float(x), float(y), float(z)))
            values.append(float(v))

    for name in names[1:]:
        with open(os.path.join(base, name)) as f:
            for line in f:
                values.append(float(line.split()[3]))

    rss = math.sqrt(sum(v * v for v in values))
    lengths.append(f"{sub:<8}{rss:16.6f}")

    if rss == 0.0:
        raise ValueError(f"{sub} vector L2 norm is zero, cannot normalize")

    with open(os.path.join(out_dir, sub), "w") as f:
        for i in range(3):
            for j in range(512):
                x, y, z = coords[j]
                v = values[i * 512 + j] / rss
                f.write(f"{x:4.1f} {y:4.1f} {z:4.1f} {v:15.12f}\n")

with open("lengths", "w") as f:
    for line in lengths:
        f.write(line + "\n")

