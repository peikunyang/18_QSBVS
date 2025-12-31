#!/usr/bin/env python3
import sys, os

pro_dir = sys.argv[1]
shot_dir = sys.argv[2]
outdir = sys.argv[3]

os.makedirs(outdir, exist_ok=True)

def read_pro(path):
    d = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            t = line.split()
            d[(t[0], t[1])] = float(t[-1])
    return d

def read_shot(path):
    d = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            t = line.split()
            d[(t[0], t[1])] = (float(t[2]), float(t[3]))
    return d

prot_order = {"comp": 0, "apo": 1}
lig_order = {"PPI-M": 0, "PPI-id": 1, "But1": 2, "But2": 3}

for tag in ["E_max_1", "E_max_10", "E_max_100"]:
    pro = read_pro(os.path.join(pro_dir, tag))

    shot_base = {
        "E_max_1": "1_max_1",
        "E_max_10": "2_max_10",
        "E_max_100": "3_max_100",
    }[tag]

    shots = {}
    for s in range(3, 8):
        shots[s] = read_shot(
            os.path.join(shot_dir, shot_base, f"energies_shots_{s}")
        )

    keys = sorted(
        pro.keys(),
        key=lambda k: (
            prot_order[k[0]],
            lig_order[k[1].rsplit("_", 2)[0]],
            int(k[1].rsplit("_", 2)[1]),
            int(k[1].rsplit("_", 2)[2]),
        )
    )

    with open(os.path.join(outdir, tag), "w") as f:
        f.write(
            f"{'':16s}"
            f"{'pro':>10s}"
            f"{'shot3':>18s}"
            f"{'shot4':>18s}"
            f"{'shot5':>18s}"
            f"{'shot6':>18s}"
            f"{'shot7':>18s}\n"
        )

        for k in keys:
            row = [
                f"{k[0]:4s}",
                f"{k[1]:10s}",
                f"{pro[k]:10.3f}",
            ]
            for s in range(3, 8):
                a, b = shots[s][k]
                row.append(f"{a:8.3f}")
                row.append(f"{b:8.3f}")
            f.write(" ".join(row) + "\n")

