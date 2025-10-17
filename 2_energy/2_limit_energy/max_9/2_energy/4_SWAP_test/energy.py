#!/usr/bin/env python3
import argparse, os, time, numpy as np, pennylane as qml
from math import log2

ap = argparse.ArgumentParser()
ap.add_argument("--rec_file", required=True)
ap.add_argument("--lig_file", required=True)
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

def load_unit(path, need, skip_norm=False):
    x = np.loadtxt(path, dtype=np.float64).reshape(-1)
    if x.size < need:
        raise ValueError(f"{path} has only {x.size} values; need >= {need}")
    x = x[:need]
    if not skip_norm:
        n = np.linalg.norm(x)
        if n == 0.0:
            raise ValueError(f"{path} first {need} values have zero norm")
        x = x / n
    return x.astype(np.float64)

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

LIG_BLOCKS = 16
u = load_unit(args.rec_file, 4096, skip_norm=args.skip_norm)
v = load_unit(args.lig_file, LIG_BLOCKS * 2048, skip_norm=args.skip_norm)

anc = 0
rec_reg = list(range(1, 13))
rec_amp = list(range(2, 13))

sel_bits = int(log2(LIG_BLOCKS))
lig_sel = list(range(13, 13 + sel_bits))
lig_amp = list(range(13 + sel_bits, 13 + sel_bits + 11))
lig_reg = lig_sel + lig_amp

max_wire = (lig_amp[-1] if lig_amp else (lig_sel[-1] if lig_sel else rec_reg[-1]))
shots = None if args.shots <= 0 else int(args.shots)
dev_kwargs = {"wires": max_wire + 1, "shots": shots}
if shots is not None and args.seed >= 0:
    dev_kwargs["seed"] = int(args.seed)

try:
    dev = qml.device("lightning.gpu", **dev_kwargs)
    print("Using device: lightning.gpu")
except Exception as e:
    print(f"GPU not available ({e}), fallback to lightning.qubit")
    dev = qml.device("lightning.qubit", **dev_kwargs)

@qml.qnode(dev, interface=None, diff_method=None)
def qnode(u, v, meas_wires):
    qml.AmplitudeEmbedding(u, wires=rec_reg, normalize=False)
    qml.AmplitudeEmbedding(v, wires=lig_reg, normalize=False)
    qml.Hadamard(anc)
    for k in range(11):
        qml.CSWAP(wires=[anc, rec_amp[k], lig_amp[k]])
    qml.Hadamard(anc)
    if shots is None:
        return qml.probs(wires=meas_wires)
    else:
        return qml.counts(wires=meas_wires)

t0 = time.perf_counter()
meas_wires = [anc, rec_reg[0]] + lig_sel
res = qnode(u, v, meas_wires)
t1 = time.perf_counter()
elapsed = t1 - t0

shape = (2, 2) + tuple(2 for _ in range(sel_bits))

if shots is None:
    probs = np.asarray(res, dtype=np.float64)
else:
    total = sum(res.values())
    size = 1
    for s in shape:
        size *= s
    probs = np.zeros(size, dtype=np.float64)
    idx_map = {}
    for i in range(size):
        bits = []
        rem = i
        for dim in shape[::-1]:
            bits.append(rem % dim)
            rem //= dim
        bits = list(reversed(bits))
        idx_map[tuple(bits)] = i
    for key, cnt in res.items():
        if isinstance(key, str):
            bits_tuple = tuple(int(b) for b in key)
        else:
            bits_tuple = tuple(int(b) for b in key)
        if len(bits_tuple) != (2 + sel_bits):
            continue
        anc_b = bits_tuple[0]; rb_b = bits_tuple[1]; sel_b = bits_tuple[2:] if sel_bits > 0 else tuple()
        full_idx = (anc_b, rb_b, *sel_b)
        i = idx_map[full_idx]
        probs[i] = cnt / total

P = probs.reshape(shape)

rec_names = rec_block_names(args.rec_order)
lig_names = lig_block_names(args.lig_order, LIG_BLOCKS)

rows = []
for rb in (0, 1):
    for lb in range(LIG_BLOCKS):
        bits = [ (lb >> j) & 1 for j in range(sel_bits-1, -1, -1) ] if sel_bits > 0 else []
        idx0 = (0, rb, *bits) if sel_bits > 0 else (0, rb)
        idx1 = (1, rb, *bits) if sel_bits > 0 else (1, rb)
        p0_joint = float(P[idx0])
        p1_joint = float(P[idx1])
        denom = p0_joint + p1_joint
        if denom <= 0.0:
            p0_cond = 0.5
            abs_ip = 0.0
        else:
            p0_cond = p0_joint / denom
            val = max(0.0, 2.0 * p0_cond - 1.0)
            abs_ip = float(np.sqrt(val))
        energy = abs_ip * args.rec_norm * args.lig_norm - args.elevate
        rows.append((rec_names[rb], lig_names[lb], p0_joint, p0_cond, abs_ip, energy, elapsed, shots if shots is not None else 0))

base_out = args.outfile
auto_out = f"{base_out}"

with open(auto_out, "w") as f:
    f.write(f"{'rec_block':30s} {'lig_block':30s} {'p0_joint':24s} {'p0_cond':24s} {'abs_ip(qm)':14s} {'位能':14s} {'elapsed_s':12s} {'shots':12s}\n")
    for rname, lname, p0j, p0c, aip, eng, el, sh in rows:
        f.write(f"{rname:30s} {lname:30s} {p0j:.20f} {p0c:.20f} {aip:14.10f} {eng:14.10f} {el:12.6f} {sh:12d}\n")

print(f"saved {len(rows)} rows to {auto_out}")

