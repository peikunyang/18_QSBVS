#!/usr/bin/env python3
# pairwise_batch.py
# Pairwise protein–ligand energies (Electrostatics with AD4 4r dielectric, and 12-6 LJ vdW).
# Batch mode over directories; strip trailing "_FF" from displayed names; aligned columns.

import argparse, os, glob, numpy as np

K_E = 332.0522173   # kcal·Å/(mol·e^2)
RMIN = 1.0e-6       # Å, avoid singularity

def read_ff(path):
    """Read FF-like table: type x y z q Rii eps per line."""
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
            q.append(  float(t[4]))
            Rii.append(float(t[5]))
            eps.append(float(t[6]))
    return (np.array(types, dtype=object),
            np.array(xyz,   dtype=float),
            np.array(q,     dtype=float),
            np.array(Rii,   dtype=float),
            np.array(eps,   dtype=float))

def pairwise_energies(pro_xyz, pro_q, pro_R, pro_eps,
                      lig_xyz, lig_q, lig_R, lig_eps):
    """Compute electrostatics (AD4 4r) and LJ 12-6 with Lorentz–Berthelot mixing."""
    if lig_xyz.size == 0 or pro_xyz.size == 0:
        return 0.0, 0.0

    dx = lig_xyz[:,None,0] - pro_xyz[None,:,0]
    dy = lig_xyz[:,None,1] - pro_xyz[None,:,1]
    dz = lig_xyz[:,None,2] - pro_xyz[None,:,2]
    r  = np.sqrt(dx*dx + dy*dy + dz*dz)
    r  = np.maximum(r, RMIN)

    qq = lig_q[:,None] * pro_q[None,:]
    E_elec = np.sum(K_E * qq / (4.0 * r * r))

    Rij   = 0.5 * (lig_R[:,None] + pro_R[None,:])
    epsij = np.sqrt(np.clip(lig_eps[:,None],0,None) * np.clip(pro_eps[None,:],0,None))
    x   = Rij / r
    x6  = x**6
    E_vdw = np.sum(4.0 * epsij * (x6*x6 - x6))

    return float(E_elec), float(E_vdw)

def collect_files(d, pattern):
    """Return sorted list of regular files in directory d matching glob pattern."""
    paths = [p for p in glob.glob(os.path.join(d, pattern)) if os.path.isfile(p)]
    return sorted(paths)

def display_stem(path):
    """Get basename without extension and drop trailing '_FF' if present."""
    base = os.path.basename(path)
    stem, _ = os.path.splitext(base)
    if stem.endswith("_FF"):
        stem = stem[:-3]
    return stem

def main():
    ap = argparse.ArgumentParser(description="Batch pairwise energies for all proteins x all ligands.")
    ap.add_argument("--protein_dir", required=True, help="Directory containing protein FF files.")
    ap.add_argument("--lig_dir",     required=True, help="Directory containing ligand FF files.")
    ap.add_argument("--protein_pattern", default="*", help="Glob for protein files (default: *)")
    ap.add_argument("--lig_pattern",     default="*", help="Glob for ligand files (default: *)")
    ap.add_argument("--out", default="energies_pairwise.txt", help="Output table file")
    args = ap.parse_args()

    pdir = os.path.abspath(args.protein_dir)
    ldir = os.path.abspath(args.lig_dir)
    if not os.path.isdir(pdir):
        raise SystemExit(f"[ERROR] Protein directory not found: {pdir}")
    if not os.path.isdir(ldir):
        raise SystemExit(f"[ERROR] Ligand directory not found: {ldir}")

    protein_files = collect_files(pdir, args.protein_pattern)
    ligand_files  = collect_files(ldir, args.lig_pattern)
    if not protein_files:
        raise SystemExit(f"[WARN] No protein files matched in: {pdir} (pattern: {args.protein_pattern})")
    if not ligand_files:
        raise SystemExit(f"[WARN] No ligand files matched in: {ldir} (pattern: {args.lig_pattern})")

    # Pre-read proteins once
    proteins = []
    protein_names = []
    for pf in protein_files:
        _, p_xyz, p_q, p_R, p_eps = read_ff(pf)
        proteins.append((pf, p_xyz, p_q, p_R, p_eps))
        protein_names.append(display_stem(pf))

    ligand_names = [display_stem(lf) for lf in ligand_files]

    # Column widths for aligned output
    col1_w = max(12, len("#protein"), max(len(n) for n in protein_names))
    col2_w = max(12, len("ligand"),  max(len(n) for n in ligand_names))
    num_w  = 18  # width for numeric columns

    with open(args.out, "w", encoding="utf-8") as w:
        w.write(f"{'#protein':<{col1_w}}  {'ligand':<{col2_w}}  {'E_elec(kcal/mol)':>{num_w}}  {'E_vdw(kcal/mol)':>{num_w}}\n")
        for (pf, p_xyz, p_q, p_R, p_eps), p_name in zip(proteins, protein_names):
            for lf, l_name in zip(ligand_files, ligand_names):
                _, l_xyz, l_q, l_R, l_eps = read_ff(lf)
                Eelec, Evdw = pairwise_energies(p_xyz, p_q, p_R, p_eps, l_xyz, l_q, l_R, l_eps)
                w.write(f"{p_name:<{col1_w}}  {l_name:<{col2_w}}  {Eelec:{num_w}.6f}  {Evdw:{num_w}.6f}\n")

if __name__ == "__main__":
    main()

