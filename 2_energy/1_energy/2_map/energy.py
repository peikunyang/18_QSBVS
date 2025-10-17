#!/usr/bin/env python3
# step3_dotmaps_batch_aligned.py
# Batch: dot receptor maps with ligand maps; aligned table with dynamic column widths.

import argparse, os, numpy as np

# ---------- I/O helpers ----------
def read_map_txt(path):
    xs, ys, zs, vals = [], [], [], []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            t = s.split()
            if len(t) < 4:
                continue
            xs.append(float(t[0])); ys.append(float(t[1])); zs.append(float(t[2])); vals.append(float(t[3]))
    xs_u = np.unique(np.array(xs))
    ys_u = np.unique(np.array(ys))
    zs_u = np.unique(np.array(zs))
    V = np.empty((len(xs_u), len(ys_u), len(zs_u)), dtype=float)
    d = {(x, y, z): v for x, y, z, v in zip(xs, ys, zs, vals)}
    for i, x in enumerate(xs_u):
        for j, y in enumerate(ys_u):
            for k, z in enumerate(zs_u):
                V[i, j, k] = d[(x, y, z)]
    return xs_u, ys_u, zs_u, V

def list_subdirs(root):
    return sorted([d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))])

def resolve_by_stem(dir_path, stem):
    """Find a file in dir whose basename (without extension) equals stem, fallback to 'contains'."""
    p0 = os.path.join(dir_path, stem)
    if os.path.isfile(p0):
        return p0
    for fn in os.listdir(dir_path):
        fp = os.path.join(dir_path, fn)
        if os.path.isfile(fp) and os.path.splitext(fn)[0] == stem:
            return fp
    for fn in os.listdir(dir_path):
        fp = os.path.join(dir_path, fn)
        if os.path.isfile(fp) and stem in os.path.splitext(fn)[0]:
            return fp
    raise FileNotFoundError(f"Cannot resolve file with stem '{stem}' in '{dir_path}'")

def load_receptor_maps(rec_dir, names=("map_elec", "map_vdw_C", "map_vdw_OA")):
    elec_path = resolve_by_stem(rec_dir, names[0])
    vdwC_path = resolve_by_stem(rec_dir, names[1])
    vdwOA_path= resolve_by_stem(rec_dir, names[2])

    xs_e, ys_e, zs_e, phi = read_map_txt(elec_path)
    xs_c, ys_c, zs_c, E_C  = read_map_txt(vdwC_path)
    xs_o, ys_o, zs_o, E_OA = read_map_txt(vdwOA_path)

    if not (np.array_equal(xs_e, xs_c) and np.array_equal(xs_e, xs_o) and
            np.array_equal(ys_e, ys_c) and np.array_equal(ys_e, ys_o) and
            np.array_equal(zs_e, zs_c) and np.array_equal(zs_e, zs_o)):
        raise ValueError(f"Receptor maps in '{rec_dir}' do not share the same grid.")
    return (xs_e, ys_e, zs_e), (phi, E_C, E_OA)

def load_ligand_maps(lig_dir):
    charge_path = resolve_by_stem(lig_dir, "map_charge")
    _, _, _, Q = read_map_txt(charge_path)

    def try_load(stem, like):
        try:
            p = resolve_by_stem(lig_dir, stem)
            return read_map_txt(p)[3]
        except FileNotFoundError:
            return np.zeros_like(like, dtype=float)

    N_C  = try_load("map_occ_C",  Q)
    N_OA = try_load("map_occ_OA", Q)
    return Q, N_C, N_OA

def format_num(x, fmt, prec):
    if fmt == "e":
        return f"{x:.{prec}e}"
    if fmt == "f":
        return f"{x:.{prec}f}"
    return f"{x:.{prec}g}"

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(description="Dot receptor maps with ligand maps over all subdirectories (aligned output).")
    ap.add_argument("--rec_map_dir", required=True, help="Root dir: each subdir has map_elec/map_vdw_C/map_vdw_OA")
    ap.add_argument("--lig_map_dir", required=True, help="Root dir: each subdir has map_charge[/map_occ_C/map_occ_OA]")
    ap.add_argument("--out", default="energies_maps.txt", help="Output table file")
    ap.add_argument("--map_elec",  default="map_elec",  help="Basename of receptor electrostatic map")
    ap.add_argument("--map_vdw_c", default="map_vdw_C", help="Basename of receptor C-type vdW map")
    ap.add_argument("--map_vdw_oa",default="map_vdw_OA",help="Basename of receptor OA-type vdW map")
    ap.add_argument("--fmt", choices=["e","f","g"], default="e", help="Numeric format: e=scientific (default), f=fixed, g=general")
    ap.add_argument("--prec", type=int, default=6, help="Precision (decimals for e/f, significant digits for g)")
    args = ap.parse_args()

    rec_root = os.path.abspath(args.rec_map_dir)
    lig_root = os.path.abspath(args.lig_map_dir)
    if not os.path.isdir(rec_root):
        raise SystemExit(f"[ERROR] Receptor map root not found: {rec_root}")
    if not os.path.isdir(lig_root):
        raise SystemExit(f"[ERROR] Ligand map root not found: {lig_root}")

    rec_names = list_subdirs(rec_root)
    lig_names = list_subdirs(lig_root)
    if not rec_names:
        raise SystemExit(f"[WARN] No receptor subfolders in: {rec_root}")
    if not lig_names:
        raise SystemExit(f"[WARN] No ligand subfolders in: {lig_root}")

    # Preload receptor maps
    rec_maps = {}
    for rec in rec_names:
        try:
            _, (phi, E_C, E_OA) = load_receptor_maps(os.path.join(rec_root, rec),
                                                     (args.map_elec, args.map_vdw_c, args.map_vdw_oa))
            rec_maps[rec] = (phi, E_C, E_OA)
        except Exception as e:
            print(f"[SKIP] Receptor '{rec}': {e}")

    # Preload ligand maps
    lig_maps = {}
    for lig in lig_names:
        try:
            lig_maps[lig] = load_ligand_maps(os.path.join(lig_root, lig))
        except Exception as e:
            print(f"[SKIP] Ligand '{lig}': {e}")

    # Compute all combinations; collect rows first (two-pass for alignment)
    rows = []
    for rec, (phi, E_C, E_OA) in rec_maps.items():
        for lig, (Q, N_C, N_OA) in lig_maps.items():
            if Q.shape != phi.shape or N_C.shape != E_C.shape or N_OA.shape != E_OA.shape:
                print(f"[SKIP] Shape mismatch: rec='{rec}', lig='{lig}' "
                      f"(Q{Q.shape} vs phi{phi.shape}, NC{N_C.shape} vs EC{E_C.shape}, NOA{N_OA.shape} vs EOA{E_OA.shape})")
                continue
            E_elec = float(np.sum(Q * phi))
            E_vdw  = float(np.sum(N_C * E_C) + np.sum(N_OA * E_OA))
            E_sum  = E_elec + E_vdw
            rows.append((rec, lig, E_elec, E_vdw, E_sum))

    # Decide column widths dynamically
    rec_w = max(12, len("#receptor"), *(len(r[0]) for r in rows)) if rows else len("#receptor")
    lig_w = max(12, len("ligand"),   *(len(r[1]) for r in rows)) if rows else len("ligand")

    # Pre-format numbers to determine widths independently for the three numeric columns
    elec_strs = [format_num(r[2], args.fmt, args.prec) for r in rows]
    vdw_strs  = [format_num(r[3], args.fmt, args.prec) for r in rows]
    sum_strs  = [format_num(r[4], args.fmt, args.prec) for r in rows]
    elec_w = max(len("E_elec(kcal/mol)"), *(len(s) for s in elec_strs)) if rows else len("E_elec(kcal/mol)")
    vdw_w  = max(len("E_vdw(kcal/mol)"),  *(len(s) for s in vdw_strs))  if rows else len("E_vdw(kcal/mol)")
    sum_w  = max(len("E_sum(kcal/mol)"),  *(len(s) for s in sum_strs))  if rows else len("E_sum(kcal/mol)")

    # Write table
    with open(args.out, "w", encoding="utf-8") as w:
        w.write(f"{'#receptor':<{rec_w}}  {'ligand':<{lig_w}}  {'E_elec(kcal/mol)':>{elec_w}}  {'E_vdw(kcal/mol)':>{vdw_w}}  {'E_sum(kcal/mol)':>{sum_w}}\n")
        for (rec, lig, E_elec, E_vdw, E_sum), se, sv, ss in zip(rows, elec_strs, vdw_strs, sum_strs):
            w.write(f"{rec:<{rec_w}}  {lig:<{lig_w}}  {se:>{elec_w}}  {sv:>{vdw_w}}  {ss:>{sum_w}}\n")

    print(f"[DONE] Wrote results to {os.path.abspath(args.out)}")

if __name__ == "__main__":
    main()

