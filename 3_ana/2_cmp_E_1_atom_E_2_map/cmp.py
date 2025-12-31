#!/usr/bin/env python3
import sys

atom_file = sys.argv[1]
map_file  = sys.argv[2]
out_file  = sys.argv[3]

W_REC = 10
W_LIG = 15
W_E10 = 12
W_E12 = 14

def fmt(x, w):
    s = "NA" if abs(x) > 9999 else f"{x:.3f}"
    return f"{s:>{w}s}"

atom = {}
with open(atom_file) as f:
    next(f)
    for line in f:
        t = line.split()
        if len(t) < 4:
            continue
        atom[(t[0], t[1])] = (float(t[2]), float(t[3]))

with open(map_file) as f, open(out_file, "w") as w:
    w.write(
        f"{'#receptor':<{W_REC}s} {'ligand':<{W_LIG}s} "
        f"{'E_elec(atom)':>{W_E10}s} {'E_elec(map)':>{W_E10}s} "
        f"{'E_vdw(atom)':>{W_E12}s} {'E_vdw(map)':>{W_E12}s} "
        f"{'E_sum(atom)':>{W_E12}s} {'E_sum(map)':>{W_E12}s}\n"
    )

    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        t = line.split()
        rec, lig_map = t[0], t[1]
        e_elec_map = float(t[2])
        e_vdw_map  = float(t[3])
        e_sum_map  = float(t[4])

        lig_atom = lig_map.replace("_", "_FF_", 1)
        if (rec, lig_atom) not in atom:
            continue

        e_elec_atom, e_vdw_atom = atom[(rec, lig_atom)]
        e_sum_atom = e_elec_atom + e_vdw_atom

        w.write(
            f"{rec:<{W_REC}s} {lig_map:<{W_LIG}s} "
            f"{fmt(e_elec_atom, W_E10)} {fmt(e_elec_map, W_E10)} "
            f"{fmt(e_vdw_atom,  W_E12)} {fmt(e_vdw_map,  W_E12)} "
            f"{fmt(e_sum_atom,  W_E12)} {fmt(e_sum_map,  W_E12)}\n"
        )

