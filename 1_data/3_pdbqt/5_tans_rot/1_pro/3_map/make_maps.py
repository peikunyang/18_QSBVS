#!/usr/bin/env python3
import argparse
import numpy as np
import os

# --------- Input: each line "type x y z q Rii eps" ---------
def read_protein(path):
    """Read receptor-like per-atom records: type x y z q Rii eps."""
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
    return (np.array(types, dtype=object),
            np.array(coords, dtype=float),
            np.array(charges, dtype=float),
            np.array(rii, dtype=float),
            np.array(eps, dtype=float))

# --------- Read ligand atom-type params (only C, OA) from AD4.1_bound.dat ---------
def load_ad4_params(ad4dat):
    """Load AD4 vdw params for ligand types; must include C and OA."""
    need = ("C", "OA")
    tab = {}
    with open(ad4dat, "r", encoding="latin-1") as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#;!":
                continue
            t = s.split()
            # e.g., "atom_par  C   4.000   0.150"
            if t[0].lower().startswith("atom_par") and len(t) >= 4:
                tab[t[1]] = (float(t[2]), float(t[3]))
    miss = [k for k in need if k not in tab]
    if miss:
        raise RuntimeError("Missing AD4 params: " + ",".join(miss))
    return tab

# --------- Grid helpers ---------
def make_grid(vmin, vmax, h):
    """Build a cubic grid from vmin..vmax with spacing h (Å)."""
    xs = np.arange(vmin, vmax + 1e-9, h)
    ys = np.arange(vmin, vmax + 1e-9, h)
    zs = np.arange(vmin, vmax + 1e-9, h)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    return xs, ys, zs, X, Y, Z

def write_map(path, xs, ys, zs, M):
    """Write scalar field M on (xs,ys,zs) grid to plain text."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as w:
        for i, x in enumerate(xs):
            for j, y in enumerate(ys):
                for k, z in enumerate(zs):
                    w.write(f"{x:8.3f} {y:8.3f} {z:8.3f} {M[i,j,k]:16.6f}\n")

# --------- Physics: AD4 4r electrostatics and LJ vdW ---------
def electrostatic_map_4r(X, Y, Z, coords, charges, k_e=332.0522173, rmin=1e-6):
    """AD4 default: epsilon(r) = 4r  =>  phi = k_e * q / (4 * r^2)."""
    phi = np.zeros_like(X, dtype=float)
    for (xi, yi, zi), qi in zip(coords, charges):
        dx = X - xi; dy = Y - yi; dz = Z - zi
        r = np.sqrt(dx*dx + dy*dy + dz*dz)
        r = np.maximum(r, rmin)
        phi += k_e * qi / (4.0 * r * r)
    return phi

def vdw_map(X, Y, Z, coords, Rrec, erec, R_lig, eps_lig, rmin=1e-6):
    """Lennard-Jones 12-6 with Lorentz-Berthelot mixing."""
    E = np.zeros_like(X, dtype=float)
    # Precompute ligand params to scalars
    Rl = float(R_lig)
    el = float(max(eps_lig, 0.0))
    for (xi, yi, zi), Rr, er in zip(coords, Rrec, erec):
        Rij   = 0.5 * (Rl + float(Rr))                     # Lorentz
        epsij = float(np.sqrt(el * max(float(er), 0.0)))   # Berthelot
        if epsij == 0.0 or Rij == 0.0:
            continue
        dx = X - xi; dy = Y - yi; dz = Z - zi
        r = np.sqrt(dx*dx + dy*dy + dz*dz)
        r = np.maximum(r, rmin)
        x = Rij / r
        x6 = x**6
        E += 4.0 * epsij * (x6*x6 - x6)
    return E

def derive_out_subdir_name(in_path):
    """Use input file stem, drop a trailing '_FF' (if present)."""
    stem = os.path.splitext(os.path.basename(in_path))[0]
    if stem.endswith("_FF"):
        stem = stem[:-3]
    return stem

def process_one_file(in_path, params, xs, ys, zs, X, Y, Z, out_root):
    """Compute maps for one input file and write to ./maps/<stem>/..."""
    types, coords, charges, Rrec, erec = read_protein(in_path)

    phi  = electrostatic_map_4r(X, Y, Z, coords, charges)       # electrostatics
    E_C  = vdw_map(X, Y, Z, coords, Rrec, erec, *params["C"])   # vdW for C
    E_OA = vdw_map(X, Y, Z, coords, Rrec, erec, *params["OA"])  # vdW for OA

    subdir = derive_out_subdir_name(in_path)
    outdir = os.path.join(out_root, subdir)

    write_map(os.path.join(outdir, "map_elec"),  xs, ys, zs, phi)
    write_map(os.path.join(outdir, "map_vdw_C"), xs, ys, zs, E_C)
    write_map(os.path.join(outdir, "map_vdw_OA"), xs, ys, zs, E_OA)
    print(f"[OK] {os.path.basename(in_path)} -> {outdir}")

def main():
    ap = argparse.ArgumentParser(
        description="Batch AD4 4r electrostatic map + vdW maps (C & OA) for all files in a directory."
    )
    # Input dir defaults to ../../1_FF/FF/1_pro relative to this script
    ap.add_argument("--indir", default=None,
                    help="Input directory (default: ../../1_FF/FF/1_pro relative to script).")
    ap.add_argument("--ad4dat", required=True,
                    help="Path to AD4.1_bound.dat (must contain atom_par for C and OA).")
    ap.add_argument("--outdir", default="maps",
                    help="Output root directory (default: ./maps).")
    ap.add_argument("--min", type=float, default=-1.2, help="Grid min (Å).")
    ap.add_argument("--max", type=float, default= 1.6, help="Grid max (Å).")
    ap.add_argument("--spacing", type=float, default=0.4, help="Grid spacing (Å).")
    args = ap.parse_args()

    # Resolve input directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if args.indir is None:
        in_dir = os.path.abspath(os.path.join(script_dir, "..", "..", "1_FF", "FF", "1_pro"))
    else:
        in_dir = os.path.abspath(args.indir)

    if not os.path.isdir(in_dir):
        raise SystemExit(f"[ERROR] Input directory not found: {in_dir}")

    # Prepare output root
    out_root = os.path.abspath(os.path.join(script_dir, args.outdir))
    os.makedirs(out_root, exist_ok=True)

    # Load AD4 ligand params (C and OA)
    params = load_ad4_params(os.path.abspath(args.ad4dat))

    # Build grid once and reuse for all inputs
    xs, ys, zs, X, Y, Z = make_grid(args.min, args.max, args.spacing)

    # Enumerate files (all regular files, ignore hidden)
    entries = sorted(os.listdir(in_dir))
    files = [os.path.join(in_dir, e) for e in entries
             if not e.startswith(".") and os.path.isfile(os.path.join(in_dir, e))]

    if not files:
        raise SystemExit(f"[WARN] No files found in: {in_dir}")

    print(f"[INFO] Input dir: {in_dir}")
    print(f"[INFO] Output root: {out_root}")
    print(f"[INFO] Files: {len(files)}")
    for fp in files:
        process_one_file(fp, params, xs, ys, zs, X, Y, Z, out_root)

    print("[DONE] All maps generated.")

if __name__ == "__main__":
    main()

