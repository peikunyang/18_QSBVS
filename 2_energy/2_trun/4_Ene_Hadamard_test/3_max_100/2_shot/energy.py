#!/usr/bin/env python3
import argparse
import os
import time
import numpy as np
import pennylane as qml
import multiprocessing as mp

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
    with open(path) as f:
        for line in f:
            if line.strip():
                k, v = line.split()
                d[k] = float(v)
    return d


def load_map_vec_2048(path):
    x = np.loadtxt(path)[:, -1]
    v = np.zeros(DIM_PAD)
    v[:DIM_IN] = x
    n = np.linalg.norm(v)
    if n == 0.0:
        raise RuntimeError(f"zero norm: {path}")
    return (v / n).astype(np.complex128)


def ctrlH(cw, cv, t):
    qml.ctrl(qml.Hadamard, control=cw, control_values=cv)(t)


def ctrlX(cw, cv, t):
    qml.ctrl(qml.PauliX, control=cw, control_values=cv)(t)


def ctrlSWAP(cw, cv, a, b):
    qml.ctrl(qml.SWAP, control=cw, control_values=cv)(wires=[a, b])


def apply_rot_ccw90_z(cw, cv):
    for a, b in zip(X_WIRES, Y_WIRES):
        ctrlSWAP(cw, cv, a, b)
    for a in X_WIRES:
        ctrlX(cw, cv, a)


def apply_shift_x_plus1(cw, cv):
    x0, x1, x2 = X_WIRES
    ctrlX(cw + [x0, x1], cv + [1, 1], x2)
    ctrlX(cw + [x0], cv + [1], x1)
    ctrlX(cw, cv, x0)


def one_run(args):
    shots, rec_vecs, lig_vecs, rec_norms, lig_norms = args

    seed = int.from_bytes(os.urandom(4), "little")
    np.random.seed(seed)

    dev = qml.device(
        "lightning.qubit",
        wires=N_WIRES,
        shots=shots,
        seed=seed
    )

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

    probs = np.asarray(qnode(*rec_vecs, *lig_vecs))
    P = probs.reshape(2, 2, 2, 2, 2, 2)

    out = {}
    for r, rec in enumerate(REC_BLOCKS):
        for l, fam in enumerate(LIG_FAMS):
            b1, b2 = (l >> 1) & 1, l & 1
            for rot in (0, 1):
                for sh in (0, 1):
                    p0 = P[0, r, b1, b2, rot, sh]
                    p1 = P[1, r, b1, b2, rot, sh]
                    ip = 0.0 if p0 + p1 == 0 else (2 * p0 / (p0 + p1) - 1)
                    en = ip * rec_norms[rec] * lig_norms[fam]
                    out[(rec, f"{fam}_{rot}_{sh}")] = en
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rec_file", required=True)
    ap.add_argument("--lig_file", required=True)
    ap.add_argument("--rec_norm", required=True)
    ap.add_argument("--lig_norm", required=True)
    ap.add_argument("--repeat", type=int, default=10)
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--nproc", type=int, default=None)
    args = ap.parse_args()

    t0 = time.time()

    rec_norms = load_norms(args.rec_norm)
    lig_norms = load_norms(args.lig_norm)

    rec_vecs = [load_map_vec_2048(os.path.join(args.rec_file, r)) for r in REC_BLOCKS]
    lig_vecs = [load_map_vec_2048(os.path.join(args.lig_file, f)) for f in LIG_FAMS]

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "energy")
    os.makedirs(outdir, exist_ok=True)

    nproc = args.nproc or min(mp.cpu_count(), args.repeat)

    for n in range(args.nmin, args.nmax + 1):
        shots = 10 ** n
        print(f"[INFO] shots = 10^{n}")

        tasks = [
            (shots, rec_vecs, lig_vecs, rec_norms, lig_norms)
            for _ in range(args.repeat)
        ]

        with mp.Pool(nproc) as pool:
            results = pool.map(one_run, tasks)

        energy_dict = {}
        for res in results:
            for k, v in res.items():
                energy_dict.setdefault(k, []).append(v)

        outfile = os.path.join(outdir, f"energies_shots_{n}")
        with open(outfile, "w") as f:
            for rec in REC_BLOCKS:
                for fam in LIG_FAMS:
                    for rot in (0, 1):
                        for sh in (0, 1):
                            key = (rec, f"{fam}_{rot}_{sh}")
                            ens = np.array(energy_dict[key])
                            mean = ens.mean()
                            std = ens.std(ddof=1) if len(ens) > 1 else 0.0

                            label = f"{fam}_{rot}_{sh}"

                            f.write(
                                f"{rec:4s} {label:10s} "
                                f"{mean:11.4e} {std:11.4e} "
                                + " ".join(f"{e:11.4e}" for e in ens)
                                + "\n"
                            )

    t1 = time.time()
    print(f"Total time: {t1 - t0:.2f} seconds")


if __name__ == "__main__":
    main()

