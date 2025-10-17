#!/usr/bin/env python3
import argparse, os, time, numpy as np, pennylane as qml
from math import log2

ap = argparse.ArgumentParser()
ap.add_argument("--rec_file", required=True, help="4096-length vector (2*2048)")
ap.add_argument("--lig_file", required=True, help="32768-length vector (16*2048)")
ap.add_argument("--outfile", required=True)
ap.add_argument("--rec_order", default=None)
ap.add_argument("--lig_order", default=None)
ap.add_argument("--rec_norm", type=float, required=True)
ap.add_argument("--lig_norm", type=float, required=True)
ap.add_argument("--shots", type=int, default=0)
ap.add_argument("--seed", type=int, default=-1)
ap.add_argument("--skip_norm", action="store_true")
ap.add_argument("--elevate", type=float, required=True)
args = ap.parse_args()

def load_vec(path, expect_len, skip_norm=False):
    x = np.loadtxt(path, dtype=np.float64).reshape(-1)
    if x.size != expect_len:
        raise ValueError(f"{path} length {x.size}, expected {expect_len}")
    x = x.astype(np.complex128)
    if not skip_norm:
        n = np.linalg.norm(x)
        if n == 0.0:
            raise ValueError(f"{path} has zero norm")
        x = x / n
    return x

u_vec = load_vec(args.rec_file, 4096, skip_norm=args.skip_norm)
v_vec = load_vec(args.lig_file, 32768, skip_norm=args.skip_norm)

anc = 0
r_idx_wires = [1]
l_idx_wires = [2,3,4,5]
data_wires  = list(range(6, 17))
all_wires = [anc] + r_idx_wires + l_idx_wires + data_wires

shots = None if args.shots <= 0 else int(args.shots)
dev_kwargs = {"wires": len(all_wires), "shots": shots}
if shots is not None and args.seed >= 0:
    dev_kwargs["seed"] = int(args.seed)

try:
    dev = qml.device("lightning.gpu", **dev_kwargs)
    print("Using device: lightning.gpu")
except Exception as e:
    print(f"GPU not available ({e}), fallback to lightning.qubit")
    dev = qml.device("lightning.qubit", **dev_kwargs)

@qml.qnode(dev, interface=None, diff_method=None)
def qnode_two_loaders(u_vec, v_vec):
    qml.Hadamard(anc)
    qml.ctrl(qml.StatePrep, control=[anc], control_values=[0])(u_vec, wires=r_idx_wires + data_wires)
    for w in l_idx_wires:
        qml.ctrl(qml.Hadamard, control=[anc], control_values=[0])(w)
    qml.ctrl(qml.StatePrep, control=[anc], control_values=[1])(v_vec, wires=l_idx_wires + data_wires)
    qml.ctrl(qml.Hadamard, control=[anc], control_values=[1])(r_idx_wires[0])
    qml.Hadamard(anc)
    meas_wires = [anc] + r_idx_wires + l_idx_wires
    if shots is None:
        return qml.probs(wires=meas_wires)
    else:
        return qml.counts(wires=meas_wires)

import time
t0 = time.perf_counter()
res = qnode_two_loaders(u_vec, v_vec)
t1 = time.perf_counter()
elapsed = t1 - t0

def bits_to_int(bits):
    v = 0
    for b in bits:
        v = (v << 1) | b
    return v

if shots is None:
    probs = np.asarray(res, dtype=np.float64)
    P = np.zeros((2,2,16), dtype=np.float64)
    for anc_b in (0,1):
        for r in (0,1):
            for l in range(16):
                l_bits = [(l>>3)&1, (l>>2)&1, (l>>1)&1, l&1]
                idx = int(f"{anc_b}{r}{l_bits[0]}{l_bits[1]}{l_bits[2]}{l_bits[3]}", 2)
                P[anc_b, r, l] = probs[idx]
else:
    total = sum(res.values())
    P = np.zeros((2,2,16), dtype=np.float64)
    for key, cnt in res.items():
        bits = tuple(int(b) for b in key) if not isinstance(key, str) else tuple(int(b) for b in key)
        anc_b = bits[0]; r_b = bits[1]; l = bits_to_int(bits[2:])
        P[anc_b, r_b, l] += cnt/total

rows = []
for r in range(2):
    for l in range(16):
        p0 = P[0, r, l]; p1 = P[1, r, l]
        denom = p0 + p1
        if denom <= 0:
            p0c = 0.5; ip = 0.0
        else:
            p0c = p0 / denom
            ip = 2.0*p0c - 1.0
        energy = ip * args.rec_norm * args.lig_norm - args.elevate
        rows.append((r, l, p0, p0c, ip, energy))

def rec_block_names(order_path):
    if order_path and os.path.isfile(order_path):
        names = [line.strip() for line in open(order_path) if line.strip()]
        if len(names) >= 2:
            return [names[0], names[1]]
    return ["rec_b0", "rec_b1"]

def lig_block_names(order_path, m):
    if order_path and os.path.isfile(order_path):
        names = [line.strip() for line in open(order_path) if line.strip()]
        if len(names) >= m:
            return names[:m]
    return [f"lig_b{i}" for i in range(m)]

rnames = rec_block_names(args.rec_order)
lnames = lig_block_names(args.lig_order, 16)

with open(args.outfile, "w") as f:
    f.write(f"# elapsed_s={elapsed:.6f}, shots={shots if shots is not None else 0}\n")
    f.write(f"{'rec_block':30s} {'lig_block':30s} {'p0_joint':24s} {'p0_cond':24s} {'ip':14s} {'位能':14s}\n")
    for r, l, p0, p0c, ip, en in rows:
        f.write(f"{rnames[r]:30s} {lnames[l]:30s} {p0:.20f} {p0c:.20f} {ip:14.10f} {en:14.10f}\n")

print("saved", len(rows), "rows to", args.outfile)

