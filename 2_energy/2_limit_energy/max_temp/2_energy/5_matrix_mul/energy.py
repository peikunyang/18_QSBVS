#!/usr/bin/env python3
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

import argparse, os, time, numpy as np, pennylane as qml

ap = argparse.ArgumentParser()
ap.add_argument("--rec_file", required=True, help="directory containing two 2048x2048 unitary files")
ap.add_argument("--lig_file", required=True, help="32768-length vector (16*2048)")
ap.add_argument("--outfile", required=True)
ap.add_argument("--rec_order", required=True)
ap.add_argument("--lig_order", required=True)
ap.add_argument("--rec_norm", type=float, required=True)
ap.add_argument("--lig_norm", type=float, required=True)
ap.add_argument("--elevate", type=float, required=True)
ap.add_argument("--shots", type=int, default=0)
ap.add_argument("--seed", type=int, default=-1)
ap.add_argument("--skip_norm", action="store_true")
args = ap.parse_args()

def load_names(path, expected_min):
    with open(path) as f:
        names = [line.strip() for line in f if line.strip()]
    if len(names) < expected_min:
        raise ValueError(f"{path}: need at least {expected_min} names, got {len(names)}")
    return names

def load_lig_vec(path, expect_len, skip_norm=False):
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

def load_two_unitaries_from_dir(dir_path, order_names):
    if not os.path.isdir(dir_path):
        raise ValueError(f"--rec_file must be a directory, got: {dir_path}")
    exts = (".pt", ".npy", ".txt", ".csv")
    files = [os.path.join(dir_path, f) for f in os.listdir(dir_path)
             if f.lower().endswith(exts) and os.path.isfile(os.path.join(dir_path, f))]
    if len(files) < 2:
        raise ValueError(f"{dir_path}: found {len(files)} candidate files, need at least 2")

    matched = []
    pool = files[:]
    for name in order_names[:2]:
        cand = [p for p in pool if name in os.path.basename(p)]
        if not cand:
            pool.sort()
            cand = [pool[0]]
        matched.append(cand[0])
        pool.remove(cand[0])

    def load_one(p):
        if p.lower().endswith(".pt"):
            import torch
            arr = torch.load(p, weights_only=True)
            arr = arr.detach().cpu().numpy() if hasattr(arr, "detach") else np.array(arr)
        elif p.lower().endswith(".npy"):
            arr = np.load(p)
        else:
            arr = np.loadtxt(p, dtype=np.float64)
        arr = np.asarray(arr, dtype=np.float64)
        if arr.shape != (2048, 2048):
            raise ValueError(f"{p}: expected shape (2048,2048), got {arr.shape}")
        return arr.astype(np.complex128)

    U0 = load_one(matched[0])
    U1 = load_one(matched[1])
    return U0, U1

rec_names = load_names(args.rec_order, 2)
lig_names = load_names(args.lig_order, 16)
U0, U1 = load_two_unitaries_from_dir(args.rec_file, rec_names)
lig_amp = load_lig_vec(args.lig_file, 32768, skip_norm=args.skip_norm)

n_qubits = 16
shots = None if args.shots == 0 else int(args.shots)
dev_kwargs = {"wires": n_qubits, "shots": shots}
if shots is not None and args.seed >= 0:
    dev_kwargs["seed"] = int(args.seed)

try:
    dev = qml.device("lightning.gpu", **dev_kwargs)
    print("Using device: lightning.gpu")
except Exception as e:
    print(f"GPU not available ({e}), fallback to lightning.qubit")
    dev = qml.device("lightning.qubit", **dev_kwargs)

lig_idx_wires = list(range(1, 5))
data_wires = list(range(5, 16))
all_wires = list(range(n_qubits))

@qml.qnode(dev, interface=None, diff_method=None)
def circuit(lig_amp, U0, U1):
    qml.Hadamard(0)
    qml.AmplitudeEmbedding(lig_amp, wires=lig_idx_wires + data_wires, normalize=True)
    qml.QubitUnitary(U0, wires=data_wires)
    V = (U1 @ U0.T).astype(np.complex128)
    qml.ctrl(qml.QubitUnitary, control=0)(V, wires=data_wires)
    if shots is None:
        return qml.probs(wires=all_wires)
    else:
        return qml.counts(wires=all_wires)

t0 = time.perf_counter()
res = circuit(lig_amp, U0, U1)
elapsed = time.perf_counter() - t0

def idx_for(c, j):
    j_bits = [(j >> 3) & 1, (j >> 2) & 1, (j >> 1) & 1, j & 1]
    bits = [c] + j_bits + [0]*11
    return int("".join(str(b) for b in bits), 2)

if shots is None:
    probs = np.asarray(res, dtype=np.float64)
    def get_prob(c, j):
        return probs[idx_for(c, j)]
else:
    total = sum(res.values())
    prob_map = {}
    for key, cnt in res.items():
        s = key if isinstance(key, str) else "".join(str(int(b)) for b in key)
        prob_map[s] = cnt / total
    def get_prob(c, j):
        s = f"{c}{j:04b}" + "0"*11
        return prob_map.get(s, 0.0)

rows = []
for r in range(2):
    for l in range(16):
        p = get_prob(r, l)
        dot_val = np.sqrt(p) if p > 0.0 else 0.0
        energy = dot_val * args.rec_norm * args.lig_norm * np.sqrt(32.0) - args.elevate
        rows.append((r, l, dot_val, energy))

rec_block_names = rec_names[:2]
lig_block_names = lig_names[:16]

with open(args.outfile, "w") as f:
    f.write(f"# elapsed_s={elapsed:.6f}, shots={args.shots}\n")
    header = (
        f"{'rec_block':25s}"
        f"{'lig_block':35s}"
        f"{'dot':>16s}"
        f"{'energy':>16s}\n"
    )
    f.write(header)
    f.write("-" * len(header) + "\n")
    for r, l, dot_val, energy in rows:
        f.write(
            f"{rec_block_names[r]:25s}"
            f"{lig_block_names[l]:35s}"
            f"{dot_val:16.10f}"
            f"{energy:16.10f}\n"
        )

print(f"saved {len(rows)} rows to {args.outfile}")

