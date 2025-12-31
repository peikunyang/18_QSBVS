#!/usr/bin/env python3
import argparse, os, glob, numpy as np

K_E = 332.0522173
RMIN = 1.0e-6

def read_ff(path):
    types=[]; xyz=[]; q=[]; Rii=[]; eps=[]
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            s=line.strip()
            if not s or s.startswith(("#",";","!")):
                continue
            t=s.split()
            if len(t) < 7:
                continue
            types.append(t[0])
            xyz.append([float(t[1]), float(t[2]), float(t[3])])
            q.append(float(t[4]))
            Rii.append(float(t[5]))
            eps.append(float(t[6]))
    return (np.array(types, dtype=object),
            np.array(xyz, dtype=float),
            np.array(q, dtype=float),
            np.array(Rii, dtype=float),
            np.array(eps, dtype=float))

def pairwise_energies(p_xyz, p_q, p_R, p_eps,
                      l_xyz, l_q, l_R, l_eps):
    if l_xyz.size == 0 or p_xyz.size == 0:
        return 0.0, 0.0

    dx = l_xyz[:,None,0] - p_xyz[None,:,0]
    dy = l_xyz[:,None,1] - p_xyz[None,:,1]
    dz = l_xyz[:,None,2] - p_xyz[None,:,2]
    r  = np.sqrt(dx*dx + dy*dy + dz*dz)
    r  = np.maximum(r, RMIN)

    qq = l_q[:,None] * p_q[None,:]
    E_elec = np.sum(K_E * qq / (4.0 * r * r))

    Rij   = 0.5 * (l_R[:,None] + p_R[None,:])
    epsij = np.sqrt(np.clip(l_eps[:,None],0,None) *
                    np.clip(p_eps[None,:],0,None))
    x  = Rij / r
    x6 = x**6
    E_vdw = np.sum(4.0 * epsij * (x6*x6 - x6))

    return float(E_elec), float(E_vdw)

def collect_files(d):
    return sorted([p for p in glob.glob(os.path.join(d, "*"))
                   if os.path.isfile(p)])

def protein_display(path):
    base = os.path.basename(path)
    if base.endswith("_FF"):
        base = base[:-3]
    return base

def ligand_group_name(path):
    base = os.path.basename(path)
    if "_FF_" in base:
        return base.split("_FF_")[0]
    return base

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein_dir", required=True)
    ap.add_argument("--lig_dir", required=True)
    ap.add_argument("--out", default="energies")
    args = ap.parse_args()

    pdir = os.path.abspath(args.protein_dir)
    ldir = os.path.abspath(args.lig_dir)

    protein_order = ["comp", "apo"]
    ligand_order  = ["PPI-M", "PPI-id", "But1", "But2"]

    protein_files = collect_files(pdir)
    ligand_files  = collect_files(ldir)

    protein_map = {protein_display(p): p for p in protein_files}
    ligand_groups = {}
    for lf in ligand_files:
        g = ligand_group_name(lf)
        ligand_groups.setdefault(g, []).append(lf)

    protein_files = [protein_map[n] for n in protein_order if n in protein_map]

    ligand_files_ordered = []
    for g in ligand_order:
        if g in ligand_groups:
            ligand_files_ordered.extend(sorted(ligand_groups[g]))

    proteins = []
    protein_names = []
    for pf in protein_files:
        _, p_xyz, p_q, p_R, p_eps = read_ff(pf)
        proteins.append((p_xyz, p_q, p_R, p_eps))
        protein_names.append(protein_display(pf))

    ligand_names = [os.path.basename(lf) for lf in ligand_files_ordered]

    col1_w = max(10, max(len(n) for n in protein_names))
    col2_w = max(18, max(len(n) for n in ligand_names))
    num_w  = 18

    with open(args.out, "w", encoding="utf-8") as w:
        w.write(f"{'protein':<{col1_w}}  {'ligand':<{col2_w}}"
                f"  {'E_elec(kcal/mol)':>{num_w}}"
                f"  {'E_vdw(kcal/mol)':>{num_w}}\n")

        for (p_xyz, p_q, p_R, p_eps), p_name in zip(proteins, protein_names):
            for lf, l_name in zip(ligand_files_ordered, ligand_names):
                _, l_xyz, l_q, l_R, l_eps = read_ff(lf)
                Eelec, Evdw = pairwise_energies(
                    p_xyz, p_q, p_R, p_eps,
                    l_xyz, l_q, l_R, l_eps
                )
                w.write(f"{p_name:<{col1_w}}  {l_name:<{col2_w}}"
                        f"  {Eelec:{num_w}.6f}"
                        f"  {Evdw:{num_w}.6f}\n")

if __name__ == "__main__":
    main()

