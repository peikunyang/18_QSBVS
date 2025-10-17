#!/usr/bin/env python3
import argparse, math, numpy as np, os
from typing import List, Tuple

# ---------- PDBQT IO (fixed columns only) ----------
def read_lines(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        return f.readlines()

def write_lines(path: str, lines: List[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)

def parse_xyz_fixed(line: str) -> Tuple[float, float, float]:
    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
    return x, y, z

def replace_xyz_fixed(line: str, xyz: Tuple[float, float, float]) -> str:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return line
    x, y, z = xyz
    l = list(line.rstrip("\n"))
    if len(l) < 54: l += [" "] * (54 - len(l))
    l[30:38] = list(f"{x:8.3f}")
    l[38:46] = list(f"{y:8.3f}")
    l[46:54] = list(f"{z:8.3f}")
    return "".join(l) + "\n"

def get_atom_indices_and_xyz(lines: List[str]) -> Tuple[List[int], np.ndarray]:
    idxs, pts = [], []
    for i, line in enumerate(lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x, y, z = parse_xyz_fixed(line)
                idxs.append(i); pts.append((x, y, z))
            except Exception:
                pass
    if not pts:
        raise RuntimeError("No ATOM/HETATM parsed from file.")
    return idxs, np.array(pts, dtype=np.float64)

def apply_transform(lines: List[str], R: np.ndarray, t: np.ndarray) -> List[str]:
    out = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                v = np.array(parse_xyz_fixed(line))
                v_new = (R @ v) + t
                out.append(replace_xyz_fixed(line, (v_new[0], v_new[1], v_new[2])))
            except Exception:
                out.append(line)
        else:
            out.append(line)
    return out

# ---------- Geometry ----------
def euler_xyz(rx, ry, rz) -> np.ndarray:
    cx, sx = math.cos(rx), math.sin(rx)
    cy, sy = math.cos(ry), math.sin(ry)
    cz, sz = math.cos(rz), math.sin(rz)
    Rx = np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]])
    Ry = np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])
    Rz = np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
    return Rz @ Ry @ Rx

def pca_basis(points: np.ndarray) -> np.ndarray:
    c = points.mean(axis=0, keepdims=True)
    A = points - c
    cov = np.cov(A.T)
    vals, vecs = np.linalg.eig(cov)
    order = np.argsort(vals)[::-1]
    vecs = vecs[:, order]
    if np.linalg.det(vecs) < 0:
        vecs[:, 2] *= -1
    return vecs

def cube_side_center(points: np.ndarray) -> Tuple[float, np.ndarray]:
    mins = points.min(axis=0); maxs = points.max(axis=0)
    side = float(np.max(maxs - mins))
    center = (mins + maxs) / 2.0
    return side, center

def search_best_rotation(points: np.ndarray, levels=3, coarse_deg=15.0, fine_factor=3.0):
    best_side = float("inf"); best_angles = (0.0, 0.0, 0.0)
    seeds = [np.eye(3), pca_basis(points).T]
    for R0 in seeds:
        side, _ = cube_side_center(points @ R0.T)
        if side < best_side:
            best_side = side; best_angles = (0.0, 0.0, 0.0)
    step = math.radians(coarse_deg)
    center_angles = list(best_angles)
    for _ in range(levels):
        rx0, ry0, rz0 = center_angles
        candidates = np.linspace(-step, step, 7)
        for drx in candidates:
            for dry in candidates:
                for drz in candidates:
                    rx, ry, rz = rx0 + drx, ry0 + dry, rz0 + drz
                    R = euler_xyz(rx, ry, rz)
                    rot = points @ R.T
                    side, _ = cube_side_center(rot)
                    if side < best_side:
                        best_side = side; best_angles = (rx, ry, rz)
        center_angles = list(best_angles); step = step / fine_factor
    R = euler_xyz(*best_angles)
    return R, best_side, best_angles

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(description="Rotate/translate ligand and protein; outputs saved to script directory, keeping input basenames.")
    ap.add_argument("--lig", required=True, help="Ligand PDBQT path (any directory)")
    ap.add_argument("--pro", required=True, help="Protein PDBQT path (any directory)")
    ap.add_argument("--levels", type=int, default=3)
    ap.add_argument("--coarse_step", type=float, default=15.0)
    ap.add_argument("--fine_factor", type=float, default=3.0)
    args = ap.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    lig_base = os.path.basename(args.lig)
    pro_base = os.path.basename(args.pro)

    # 若兩者檔名相同，避免覆蓋：自動加前綴
    if lig_base == pro_base:
        out_lig = os.path.join(script_dir, f"lig_{lig_base}")
        out_pro = os.path.join(script_dir, f"pro_{pro_base}")
    else:
        out_lig = os.path.join(script_dir, lig_base)
        out_pro = os.path.join(script_dir, pro_base)

    lig_lines = read_lines(args.lig)
    pro_lines = read_lines(args.pro)

    _, lig_pts = get_atom_indices_and_xyz(lig_lines)

    R, _, ang = search_best_rotation(
        lig_pts, levels=args.levels, coarse_deg=args.coarse_step, fine_factor=args.fine_factor
    )
    rot = lig_pts @ R.T
    _, center1 = cube_side_center(rot)
    t = -center1

    lig_lines_1 = apply_transform(lig_lines, R, t)
    pro_lines_1 = apply_transform(pro_lines, R, t)

    _, lig_pts_after = get_atom_indices_and_xyz(lig_lines_1)
    _, center2 = cube_side_center(lig_pts_after)
    t_resid = -center2
    if np.linalg.norm(t_resid) > 1e-6:
        lig_lines_1 = apply_transform(lig_lines_1, np.eye(3), t_resid)
        pro_lines_1 = apply_transform(pro_lines_1, np.eye(3), t_resid)
        _, lig_pts_after = get_atom_indices_and_xyz(lig_lines_1)

    write_lines(out_lig, lig_lines_1)
    write_lines(out_pro, pro_lines_1)

    mins = lig_pts_after.min(axis=0); maxs = lig_pts_after.max(axis=0)
    print("=== Summary ===")
    print(f"Saved: {out_lig}")
    print(f"Saved: {out_pro}")
    print(f"Best Euler angles (deg): rx={math.degrees(ang[0]):.3f}, ry={math.degrees(ang[1]):.3f}, rz={math.degrees(ang[2]):.3f}")
    print(f"Ligand min xyz: ({mins[0]:.3f}, {mins[1]:.3f}, {mins[2]:.3f})")
    print(f"Ligand max xyz: ({maxs[0]:.3f}, {maxs[1]:.3f}, {maxs[2]:.3f})")

if __name__ == "__main__":
    main()

