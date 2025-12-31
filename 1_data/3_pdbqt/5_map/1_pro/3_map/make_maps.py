#!/usr/bin/env python3
import argparse
import numpy as np
import os

# --------- Input: each line "type x y z q Rii eps" ---------
def read_protein(path):
    types, coords, charges, rii, eps = [], [], [], [], []
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#;!":
                continue
            t = s.split()
            if len(t) < 7:
                continue
            types.append(t[0])
            coords.append([float(t[1]), float(t[2]), float(t[3])])
            charges.append(float(t[4]))
            rii.append(float(t[5]))
            eps.append(float(t[6]))
    return (
        np.array(types, dtype=object),
        np.array(coords, dtype=float),
        np.array(charges, dtype=float),
        np.array(rii, dtype=float),
        np.array(eps, dtype=float),
    )

# --------- Read ligand atom-type params (only C, OA) ---------
def load_ad4_params(ad4dat):
    need = ("C", "OA")
    tab = {}
    with open(ad4dat, "r", encoding="latin-1") as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#;!":
                continue
            t = s.split()
            if t[0].lower().startswith("atom_par") and len(t) >= 4:
                tab[t[1]] = (float(t[2]), float(t[3]))
    miss = [k for k in need if k not in tab]
    if miss:
        raise RuntimeError("Missing AD4 params: " + ",".join(miss))
    return tab

# --------- Grid helpers ---------
def make_grid(vmin, vmax, h):
    xs = np.arange(vmin, vmax + 1e-9, h)
    ys = np.arange(vmin, vmax + 1e-9, h)
    zs = np.arange(vmin, vmax + 1e-9, h)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    return xs, ys, zs, X, Y, Z

# --------- WRITE MAP (x fastest, output z y x) ---------
def write_map(path, xs, ys, zs, M):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as w:
        for k, z in enumerate(zs):          # z slowest
            for j, y in enumerate(ys):      # y middle
                for i, x in enumerate(xs):  # x fastest
                    w.write(
                        f"{z:8.3f} {y:8.3f} {x:8.3f} {M[i,j,k]:16.6f}\n"
                    )

# --------- Physics: AD4 4r electrostatics ---------
def electrostatic_map_4r(X, Y, Z, coords, charges, k_e=332.0522173, rmin=1e-6):
    phi = np.zeros_like(X, dtype=float)
    for (xi, yi, zi), qi in zip(coords, charges):
        dx = X - xi
        dy = Y - yi
        dz = Z - zi
        r = np.sqrt(dx*dx + dy*dy + dz*dz)
        r = np.maximum(r, rmin)
        phi += k_e * qi / (4.0 * r * r)
    return phi

# --------- Physics: vdW LJ 12-6 ---------
def vdw_map(X, Y, Z, coords, Rrec, erec, R_lig, eps_lig, rmin=1e-6):
    E = np.zeros_like(X, dtype=float)
    Rl = float(R_lig)
    el = float(max(eps_lig, 0.0))
    for (xi, yi, zi), Rr, er in zip(coords, Rrec, erec):
        Rij = 0.5 * (Rl + float(Rr))
        epsij = np.sqrt(el * max(float(er), 0.0))
        if epsij == 0.0 or Rij == 0.0:
            continue
        dx = X - xi
        dy = Y - yi
        dz = Z - zi
        r = np.sqrt(dx*dx + dy*dy + dz*dz)
        r = np.maximum(r, rmin)
        x = Rij / r
        x6 = x**6
        E += 4.0 * epsij * (x6*x6 - x6)
    return E

def derive_out_subdir_name(in_path):
    stem = os.path.splitext(os.path.basename(in_path))[0]
    if stem.endswith("_FF"):
        stem = stem[:-3]
    return stem

def process_one_file(in_path, params, xs, ys, zs, X, Y, Z, out_root):
    types, coords, charges, Rrec, erec = read_protein(in_path)

    phi  = electrostatic_map_4r(X, Y, Z, coords, charges)
    E_C  = vdw_map(X, Y, Z, coords, Rrec, erec, *params["C"])
    E_OA = vdw_map(X, Y, Z, coords, Rrec, erec, *params["OA"])

    subdir = derive_out_subdir_name(in_path)
    outdir = os.path.join(out_root, subdir)

    write_map(os.path.join(outdir, "map_elec"),  xs, ys, zs, phi)
    write_map(os.path.join(outdir, "map_vdw_C"), xs, ys, zs, E_C)
    write_map(os.path.join(outdir, "map_vdw_OA"), xs, ys, zs, E_OA)

    print(f"[OK] {os.path.basename(in_path)} -> {outdir}")

def main():
    ap = argparse.ArgumentParser(
        description="Batch AD4 4r electrostatic + vdW maps (C/OA), x-fastest, output z y x."
    )
    ap.add_argument("--indir", default=None)
    ap.add_argument("--ad4dat", required=True)
    ap.add_argument("--outdir", default="maps")
    ap.add_argument("--min", type=float, default=-1.2)
    ap.add_argument("--max", type=float, default=1.6)
    ap.add_argument("--spacing", type=float, default=0.4)
    args = ap.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    if args.indir is None:
        in_dir = os.path.abspath(os.path.join(script_dir, "..", "..", "1_FF", "FF", "1_pro"))
    else:
        in_dir = os.path.abspath(args.indir)

    if not os.path.isdir(in_dir):
        raise SystemExit(f"[ERROR] Input directory not found: {in_dir}")

    out_root = os.path.abspath(os.path.join(script_dir, args.outdir))
    os.makedirs(out_root, exist_ok=True)

    params = load_ad4_params(os.path.abspath(args.ad4dat))
    xs, ys, zs, X, Y, Z = make_grid(args.min, args.max, args.spacing)

    files = [
        os.path.join(in_dir, f)
        for f in sorted(os.listdir(in_dir))
        if not f.startswith(".") and os.path.isfile(os.path.join(in_dir, f))
    ]

    if not files:
        raise SystemExit("[WARN] No input files found.")

    for fp in files:
        process_one_file(fp, params, xs, ys, zs, X, Y, Z, out_root)

    print("[DONE] All maps generated.")

if __name__ == "__main__":
    main()

