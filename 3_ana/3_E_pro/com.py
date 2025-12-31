#!/usr/bin/env python3
import sys
import math

e1_atom = sys.argv[1]
e2_map = sys.argv[2]
emax_dir = sys.argv[3]
had_dir = sys.argv[4]
outfile = sys.argv[5]

def fmt(v, width):
    if v is None:
        return f"{'NA':>{width}s}"
    if abs(v) > 999.99 or math.isnan(v):
        return f"{'NA':>{width}s}"
    return f"{v:{width}.6f}"

atom_sum = {}
with open(e1_atom) as f:
    next(f)
    for line in f:
        if not line.strip():
            continue
        t = line.split()
        lig = t[1].replace("_FF_", "_").replace("_FF", "")
        key = (t[0], lig)
        atom_sum[key] = float(t[2]) + float(t[3])

map_sum = {}
with open(e2_map) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        t = line.split()
        key = (t[0], t[1])
        map_sum[key] = float(t[4])

emax = {}
for tag in ["E_max_1", "E_max_10", "E_max_100"]:
    with open(f"{emax_dir}/{tag}") as f:
        for line in f:
            if not line.strip():
                continue
            t = line.split()
            key = (t[0], t[1])
            emax.setdefault(key, {})[tag] = float(t[2])

had = {}
for tag in ["E_max_1", "E_max_10", "E_max_100"]:
    with open(f"{had_dir}/{tag}") as f:
        for line in f:
            if not line.strip():
                continue
            t = line.split()
            try:
                val = float(t[-1])
            except ValueError:
                continue
            key = (t[0], t[1])
            had.setdefault(key, {})[tag] = val

keys = set()
keys |= atom_sum.keys()
keys |= map_sum.keys()
keys |= emax.keys()
keys |= had.keys()

prot_order = {"comp": 0, "apo": 1}
lig_order = {"PPI-M": 0, "PPI-id": 1, "But1": 2, "But2": 3}

keys = sorted(
    keys,
    key=lambda k: (
        prot_order.get(k[0], 99),
        lig_order.get(k[1].rsplit("_", 2)[0], 99),
        int(k[1].rsplit("_", 2)[1]),
        int(k[1].rsplit("_", 2)[2]),
    )
)

with open(outfile, "w") as f:
    f.write(
        " " * 18
        + f"{'atom':>12s} "
        + f"{'map':>12s} "
        + f"{'max_1':>12s} "
        + f"{'max_10':>12s} "
        + f"{'max_100':>12s} "
        + f"{'had_1':>12s} "
        + f"{'had_10':>12s} "
        + f"{'had_100':>12s}\n"
    )

    for k in keys:
        f.write(
            f"{k[0]:4s} {k[1]:12s} "
            f"{fmt(atom_sum.get(k), 12)} "
            f"{fmt(map_sum.get(k), 12)} "
            f"{fmt(emax.get(k, {}).get('E_max_1'), 12)} "
            f"{fmt(emax.get(k, {}).get('E_max_10'), 12)} "
            f"{fmt(emax.get(k, {}).get('E_max_100'), 12)} "
            f"{fmt(had.get(k, {}).get('E_max_1'), 12)} "
            f"{fmt(had.get(k, {}).get('E_max_10'), 12)} "
            f"{fmt(had.get(k, {}).get('E_max_100'), 12)}\n"
        )

