#!/usr/bin/env python3
import argparse

def key_fields(line):
    return (
        line[12:16],
        line[17:20],
        line[21:22],
        line[22:26],
    )

def parse_xyz(line):
    return float(line[30:38]), float(line[38:46]), float(line[46:54])

def replace_xyz(line, xyz):
    x, y, z = xyz
    l = list(line.rstrip("\n"))
    if len(l) < 54:
        l += [" "] * (54 - len(l))
    l[30:38] = list(f"{x:8.3f}")
    l[38:46] = list(f"{y:8.3f}")
    l[46:54] = list(f"{z:8.3f}")
    return "".join(l) + "\n"

ap = argparse.ArgumentParser()
ap.add_argument("--orig", required=True)
ap.add_argument("--chim", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

with open(args.orig) as f:
    orig = f.readlines()

with open(args.chim) as f:
    chim = f.readlines()

oi = 0
out = []

for cl in chim:
    if cl.startswith(("ATOM", "HETATM")):
        while oi < len(orig) and not orig[oi].startswith(("ATOM", "HETATM")):
            out.append(orig[oi])
            oi += 1
        if oi >= len(orig):
            raise RuntimeError("Original file shorter than Chimera file")
        if key_fields(cl) != key_fields(orig[oi]):
            raise RuntimeError(
                f"Mismatch at atom:\nCHIMERA: {cl}\nORIGINAL: {orig[oi]}"
            )
        xyz = parse_xyz(cl)
        out.append(replace_xyz(orig[oi], xyz))
        oi += 1
    else:
        out.append(cl)

with open(args.out, "w") as f:
    f.writelines(out)

print("OK: pdbqt restored with verified atom/residue consistency")

