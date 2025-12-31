#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pennylane as qml

ap = argparse.ArgumentParser()
ap.add_argument("--rec_file", required=True)
ap.add_argument("--lig_file", required=True)
ap.add_argument("--rec_norm", required=True)
ap.add_argument("--lig_norm", required=True)
ap.add_argument("--outfile", required=True)
args = ap.parse_args()

REC_BLOCKS = ["comp", "apo"]
LIG_FAMS = ["PPI-M", "PPI-id", "But1", "But2"]

DIM_IN = 1536
DIM_PAD = 2048

REC_WIRE = 0
LIG_WIRES = [1, 2]
ROT_WIRE = 3
SHIFT_WIRE = 4
DATA_WIRES = list(range(5, 16))
ANC = 16
N_WIRES = 17

Z_WIRES = [9, 8, 7]
Y_WIRES = [12, 11, 10]
X_WIRES = [15, 14, 13]

def load_norms(path):
    d = {}
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            k, v = line.split()
            d[k] = float(v)
    return d

def load_map_vec_2048(path):
    x = np.loadtxt(path, dtype=np.float64)[:, -1]
    if x.size != DIM_IN:
        raise ValueError(f"{path}: expect {DIM_IN}, got {x.size}")
    v = np.zeros(DIM_PAD, dtype=np.float64)
    v[:DIM_IN] = x
    n = np.linalg.norm(v)
    if n == 0.0:
        raise ValueError(f"{path}: zero norm")
    return (v / n).astype(np.complex128)

rec_norms = load_norms(args.rec_norm)
lig_norms = load_norms(args.lig_norm)

rec_vecs = []
for rec in REC_BLOCKS:
    rec_vecs.append(load_map_vec_2048(os.path.join(args.rec_file, rec)))

lig_vecs = []
for fam in LIG_FAMS:
    lig_vecs.append(load_map_vec_2048(os.path.join(args.lig_file, fam)))

dev = qml.device("lightning.qubit", wires=N_WIRES, shots=None)

def ctrlH(ctrl_wires, ctrl_vals, target):
    qml.ctrl(qml.Hadamard, control=ctrl_wires, control_values=ctrl_vals)(target)

def ctrlX(ctrl_wires, ctrl_vals, target):
    qml.ctrl(qml.PauliX, control=ctrl_wires, control_values=ctrl_vals)(target)

def ctrlSWAP(ctrl_wires, ctrl_vals, a, b):
    qml.ctrl(qml.SWAP, control=ctrl_wires, control_values=ctrl_vals)(wires=[a, b])

def apply_rot_ccw90_z(ctrl_wires, ctrl_vals):
    for a, b in zip(X_WIRES, Y_WIRES):
        ctrlSWAP(ctrl_wires, ctrl_vals, a, b)
    for a in X_WIRES:
        ctrlX(ctrl_wires, ctrl_vals, a)

def apply_shift_x_plus1(ctrl_wires, ctrl_vals):
    c = list(ctrl_wires)
    v = list(ctrl_vals)
    x0, x1, x2 = X_WIRES
    ctrlX(c + [x0, x1], v + [1, 1], x2)
    ctrlX(c + [x0], v + [1], x1)
    ctrlX(c, v, x0)

@qml.qnode(dev, interface=None, diff_method=None)
def qnode(rec0, rec1, lig0, lig1, lig2, lig3):
    qml.Hadamard(ANC)

    ctrlH([ANC], [0], REC_WIRE)
    for w in LIG_WIRES:
        ctrlH([ANC], [0], w)
    ctrlH([ANC], [0], ROT_WIRE)
    ctrlH([ANC], [0], SHIFT_WIRE)

    qml.ctrl(qml.StatePrep, control=[ANC, REC_WIRE], control_values=[0, 0])(rec0, wires=DATA_WIRES)
    qml.ctrl(qml.StatePrep, control=[ANC, REC_WIRE], control_values=[0, 1])(rec1, wires=DATA_WIRES)

    ctrlH([ANC], [1], REC_WIRE)
    for w in LIG_WIRES:
        ctrlH([ANC], [1], w)
    ctrlH([ANC], [1], ROT_WIRE)
    ctrlH([ANC], [1], SHIFT_WIRE)

    qml.ctrl(qml.StatePrep, control=[ANC, LIG_WIRES[0], LIG_WIRES[1]], control_values=[1, 0, 0])(lig0, wires=DATA_WIRES)
    qml.ctrl(qml.StatePrep, control=[ANC, LIG_WIRES[0], LIG_WIRES[1]], control_values=[1, 0, 1])(lig1, wires=DATA_WIRES)
    qml.ctrl(qml.StatePrep, control=[ANC, LIG_WIRES[0], LIG_WIRES[1]], control_values=[1, 1, 0])(lig2, wires=DATA_WIRES)
    qml.ctrl(qml.StatePrep, control=[ANC, LIG_WIRES[0], LIG_WIRES[1]], control_values=[1, 1, 1])(lig3, wires=DATA_WIRES)

    apply_rot_ccw90_z([ANC, ROT_WIRE], [1, 1])
    apply_shift_x_plus1([ANC, SHIFT_WIRE], [1, 1])

    qml.Hadamard(ANC)

    return qml.probs(wires=[ANC, REC_WIRE, LIG_WIRES[0], LIG_WIRES[1], ROT_WIRE, SHIFT_WIRE])

probs = np.asarray(
    qnode(rec_vecs[0], rec_vecs[1], lig_vecs[0], lig_vecs[1], lig_vecs[2], lig_vecs[3]),
    dtype=np.float64
)

P = np.zeros((2, 2, 4, 2, 2), dtype=np.float64)
for anc_b in (0, 1):
    for r in (0, 1):
        for l in range(4):
            for rot in (0, 1):
                for sh in (0, 1):
                    b1 = (l >> 1) & 1
                    b2 = l & 1
                    idx = int(f"{anc_b}{r}{b1}{b2}{rot}{sh}", 2)
                    P[anc_b, r, l, rot, sh] = probs[idx]

rows = []
for r in range(2):
    rec_name = REC_BLOCKS[r]
    rnorm = rec_norms[rec_name]
    for l in range(4):
        fam = LIG_FAMS[l]
        lnorm = lig_norms[fam]
        for rot in (0, 1):
            for sh in (0, 1):
                p0 = P[0, r, l, rot, sh]
                p1 = P[1, r, l, rot, sh]
                denom = p0 + p1
                if denom <= 0.0:
                    p0c = 0.5
                    ip = 0.0
                else:
                    p0c = p0 / denom
                    ip = 2.0 * p0c - 1.0
                en = ip * rnorm * lnorm
                liglab = f"{fam}_{rot}_{sh}"
                rows.append((rec_name, liglab, p0c, ip, en))

rec_rank = {k: i for i, k in enumerate(REC_BLOCKS)}
lig_rank = {k: i for i, k in enumerate(LIG_FAMS)}
rows.sort(key=lambda t: (rec_rank[t[0]], lig_rank[t[1].split("_")[0]], int(t[1].split("_")[1]), int(t[1].split("_")[2])))

with open(args.outfile, "w") as f:
    for rec, liglab, p0c, ip, en in rows:
        f.write(f"{rec:4s} {liglab:10s} {p0c:11.4e} {ip:11.4e} {en:11.4e}\n")

print(f"saved {len(rows)} rows to {args.outfile}")

