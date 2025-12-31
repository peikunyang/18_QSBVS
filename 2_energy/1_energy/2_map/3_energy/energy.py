#!/usr/bin/env python3
import argparse, os, numpy as np

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
            xs.append(float(t[0]))
            ys.append(float(t[1]))
            zs.append(float(t[2]))
            vals.append(float(t[3]))
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
    raise FileNotFoundError

def load_receptor_maps(rec_dir, names):
    elec_path = resolve_by_stem(rec_dir, names[0])
    vdwC_path = resolve_by_stem(rec_dir, names[1])
    vdwOA_path = resolve_by_stem(rec_dir, names[2])
    xs_e, ys_e, zs_e, phi = read_map_txt(elec_path)
    xs_c, ys_c, zs_c, E_C = read_map_txt(vdwC_path)
    xs_o, ys_o, zs_o, E_OA = read_map_txt(vdwOA_path)
    if not (np.array_equal(xs_e, xs_c) and np.array_equal(xs_e, xs_o) and
            np.array_equal(ys_e, ys_c) and np.array_equal(ys_e, ys_o) and
            np.array_equal(zs_e, zs_c) and np.array_equal(zs_e, zs_o)):
        raise ValueError
    return phi, E_C, E_OA

def load_ligand_maps(lig_dir):
    charge_path = resolve_by_stem(lig_dir, "map_charge")
    _, _, _, Q = read_map_txt(charge_path)
    def try_load(stem, like):
        try:
            return read_map_txt(resolve_by_stem(lig_dir, stem))[3]
        except:
            return np.zeros_like(like)
    N_C = try_load("map_occ_C", Q)
    N_OA = try_load("map_occ_OA", Q)
    return Q, N_C, N_OA

def format_num(x, fmt, prec):
    if fmt == "e":
        return f"{x:.{prec}e}"
    if fmt == "f":
        return f"{x:.{prec}f}"
    return f"{x:.{prec}g}"

def ligand_main_name(lig):
    return lig.split("_", 1)[0]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rec_map_dir", required=True)
    ap.add_argument("--lig_map_dir", required=True)
    ap.add_argument("--out", default="energies")
    ap.add_argument("--map_elec", default="map_elec")
    ap.add_argument("--map_vdw_c", default="map_vdw_C")
    ap.add_argument("--map_vdw_oa", default="map_vdw_OA")
    ap.add_argument("--fmt", choices=["e","f","g"], default="e")
    ap.add_argument("--prec", type=int, default=6)
    args = ap.parse_args()

    rec_root = os.path.abspath(args.rec_map_dir)
    lig_root = os.path.abspath(args.lig_map_dir)

    rec_names = list_subdirs(rec_root)
    lig_names = list_subdirs(lig_root)

    rec_maps = {}
    for rec in rec_names:
        try:
            rec_maps[rec] = load_receptor_maps(
                os.path.join(rec_root, rec),
                (args.map_elec, args.map_vdw_c, args.map_vdw_oa)
            )
        except:
            pass

    lig_maps = {}
    for lig in lig_names:
        try:
            lig_maps[lig] = load_ligand_maps(os.path.join(lig_root, lig))
        except:
            pass

    rows = []
    for rec, (phi, E_C, E_OA) in rec_maps.items():
        for lig, (Q, N_C, N_OA) in lig_maps.items():
            if Q.shape != phi.shape:
                continue
            E_elec = float(np.sum(Q * phi))
            E_vdw = float(np.sum(N_C * E_C) + np.sum(N_OA * E_OA))
            rows.append((rec, lig, E_elec, E_vdw, E_elec + E_vdw))

    protein_order = ["comp", "apo"]
    ligand_order = ["PPI-M", "PPI-id", "But1", "But2"]
    protein_rank = {p: i for i, p in enumerate(protein_order)}
    ligand_rank = {l: i for i, l in enumerate(ligand_order)}

    rows.sort(
        key=lambda r: (
            protein_rank.get(r[0], 999),
            ligand_rank.get(ligand_main_name(r[1]), 999),
            r[1]
        )
    )

    rec_w = max(len("#receptor"), *(len(r[0]) for r in rows))
    lig_w = max(len("ligand"), *(len(r[1]) for r in rows))

    elec_s = [format_num(r[2], args.fmt, args.prec) for r in rows]
    vdw_s  = [format_num(r[3], args.fmt, args.prec) for r in rows]
    sum_s  = [format_num(r[4], args.fmt, args.prec) for r in rows]

    elec_w = max(len("E_elec"), *(len(s) for s in elec_s))
    vdw_w  = max(len("E_vdw"),  *(len(s) for s in vdw_s))
    sum_w  = max(len("E_sum"),  *(len(s) for s in sum_s))

    with open(args.out, "w", encoding="utf-8") as w:
        w.write(f"{'#receptor':<{rec_w}}  {'ligand':<{lig_w}}  {'E_elec':>{elec_w}}  {'E_vdw':>{vdw_w}}  {'E_sum':>{sum_w}}\n")
        for (rec, lig, _, _, _), se, sv, ss in zip(rows, elec_s, vdw_s, sum_s):
            w.write(f"{rec:<{rec_w}}  {lig:<{lig_w}}  {se:>{elec_w}}  {sv:>{vdw_w}}  {ss:>{sum_w}}\n")

if __name__ == "__main__":
    main()

