import os
import argparse
import shutil
import sys

def rotate_90ccw(x, y, z):
    return -y, x, z

def read_map(path):
    xs, ys, zs, os_ = [], [], [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            t = line.split()
            if len(t) < 4:
                continue
            xs.append(float(t[0]))
            ys.append(float(t[1]))
            zs.append(float(t[2]))
            os_.append(float(t[3]))
    return xs, ys, zs, os_

def write_map(path, xs, ys, zs, os_):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as w:
        for x, y, z, o in zip(xs, ys, zs, os_):
            w.write(f"{x:8.3f} {y:8.3f} {z:8.3f} {o:14.6f}\n")

def idx_from_coord(val, xmin, step):
    step_m = int(round(step * 1000))
    xmin_m = int(round(xmin * 1000))
    val_m  = int(round(val  * 1000))
    k = (val_m - xmin_m + step_m // 2) // step_m
    return int(k)

def scatter_to_grid(xs, ys, zs, os_, xmin, xmax, step, rotate=False):
    n = int(round((xmax - xmin) / step)) + 1
    grid = [0.0] * (n * n * n)
    def idx3(ix, iy, iz):
        return ix * (n * n) + iy * n + iz
    cnt_in = 0
    cnt_put = 0
    for x, y, z, o in zip(xs, ys, zs, os_):
        xr, yr, zr = (-y, x, z) if rotate else (x, y, z)
        ix = idx_from_coord(xr, xmin, step)
        iy = idx_from_coord(yr, xmin, step)
        iz = idx_from_coord(zr, xmin, step)
        cnt_in += 1
        if 0 <= ix < n and 0 <= iy < n and 0 <= iz < n:
            p = idx3(ix, iy, iz)
            grid[p] = o
            cnt_put += 1
    return grid, n, cnt_in, cnt_put

def sample_from_grid(xs, ys, zs, grid, xmin, xmax, step, n):
    def idx3(ix, iy, iz):
        return ix * (n * n) + iy * n + iz
    out = []
    for x, y, z in zip(xs, ys, zs):
        ix = idx_from_coord(x, xmin, step)
        iy = idx_from_coord(y, xmin, step)
        iz = idx_from_coord(z, xmin, step)
        if 0 <= ix < n and 0 <= iy < n and 0 <= iz < n:
            out.append(grid[idx3(ix, iy, iz)])
        else:
            out.append(0.0)
    return out

def self_check_rot0(xs, ys, zs, os_, xmin, xmax, step, tol=1e-12, max_report=5):
    grid0, n, _, _ = scatter_to_grid(xs, ys, zs, os_, xmin, xmax, step, rotate=False)
    os_chk = sample_from_grid(xs, ys, zs, grid0, xmin, xmax, step, n)
    bad = []
    for i, (a, b) in enumerate(zip(os_, os_chk)):
        if abs(a - b) > tol:
            bad.append((i, a, b, xs[i], ys[i], zs[i]))
            if len(bad) >= max_report:
                break
    return bad

parser = argparse.ArgumentParser()
parser.add_argument("--input_root", required=True)
parser.add_argument("--output_root", required=True)
parser.add_argument("--xmin", type=float, default=-1.8)
parser.add_argument("--xmax", type=float, default=2.4)
parser.add_argument("--step", type=float, default=0.6)
args = parser.parse_args()

in_root = os.path.abspath(args.input_root)
out_root = os.path.abspath(args.output_root)

files = []
for root, _, fnames in os.walk(in_root):
    for fn in fnames:
        p = os.path.join(root, fn)
        if os.path.isfile(p):
            rel = os.path.relpath(p, in_root)
            top = rel.split(os.sep, 1)[0]
            sub_rel = "" if top == rel else rel[len(top)+1:]
            files.append((p, top, sub_rel))
files.sort()

tops = sorted({top for _, top, _ in files})
for top in tops:
    os.makedirs(os.path.join(out_root, top + "_rot0"), exist_ok=True)
    os.makedirs(os.path.join(out_root, top + "_rot1"), exist_ok=True)

for in_path, top, sub_rel in files:
    out_path0 = os.path.join(out_root, top + "_rot0", sub_rel) if sub_rel else os.path.join(out_root, top + "_rot0", os.path.basename(in_path))
    os.makedirs(os.path.dirname(out_path0), exist_ok=True)
    shutil.copy2(in_path, out_path0)

    xs, ys, zs, os_ = read_map(in_path)
    diffs = self_check_rot0(xs, ys, zs, os_, args.xmin, args.xmax, args.step)
    if diffs:
        print("[ERROR] rot0 self-check failed:", in_path)
        for i, a, b, x, y, z in diffs:
            print(f"  line={i+1} coord=({x:.3f},{y:.3f},{z:.3f}) orig={a:.12f} got={b:.12f}")
        sys.exit(1)

    grid1, n, c_in1, c_put1 = scatter_to_grid(xs, ys, zs, os_, args.xmin, args.xmax, args.step, rotate=True)
    os_out1 = sample_from_grid(xs, ys, zs, grid1, args.xmin, args.xmax, args.step, n)
    out_path1 = os.path.join(out_root, top + "_rot1", sub_rel) if sub_rel else os.path.join(out_root, top + "_rot1", os.path.basename(in_path))
    write_map(out_path1, xs, ys, zs, os_out1)

    print(f"[OK] {in_path} -> rot0(copy), rot1(90° kept:{c_put1}/{c_in1})")

