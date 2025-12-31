#!/usr/bin/env python3
import os

def rot(x, y, z, m):
    if m == 0:
        return x, y, z
    return -y, x, z

in_dir = "../1_ori"
out_dir = "pdbqt"
os.makedirs(out_dir, exist_ok=True)

for fn in os.listdir(in_dir):
    in_path = os.path.join(in_dir, fn)
    if not os.path.isfile(in_path):
        continue

    with open(in_path, "r") as f:
        lines = f.readlines()

    for m in (0, 90):
        out_path = os.path.join(out_dir, fn + ("_0" if m == 0 else "_1"))
        with open(out_path, "w") as fw:
            for line in lines:
                if not line.strip():
                    fw.write(line)
                    continue

                p = line.split()
                if len(p) < 7:
                    fw.write(line)
                    continue

                typ = p[0]
                x = float(p[1])
                y = float(p[2])
                z = float(p[3])
                q = float(p[4])
                r_str = p[5]
                e_str = p[6]

                xr, yr, zr = rot(x, y, z, m)

                fw.write(
                    f"{typ:>3} {xr:12.3f} {yr:12.3f} {zr:12.3f} {q:9.3f} {r_str:>7} {e_str:>12}\n"
                )

