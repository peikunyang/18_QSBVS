import os
import argparse

def read_map(path):
    xs, ys, zs, os_ = [], [], [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            t = line.split()
            if len(t) >= 4:
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

def idx3(ix, iy, iz, n):
    return ix * (n * n) + iy * n + iz

def translate_and_fill(xs, ys, zs, os_, dx, xmin, xmax, step):
    n = int(round((xmax - xmin) / step)) + 1
    grid = [0.0] * (n * n * n)
    for x, y, z, o in zip(xs, ys, zs, os_):
        xr = x + dx
        if xmin <= xr <= xmax:
            ix = int(round((xr - xmin) / step))
            iy = int(round((y - xmin) / step))
            iz = int(round((z - xmin) / step))
            if 0 <= ix < n and 0 <= iy < n and 0 <= iz < n:
                grid[idx3(ix, iy, iz, n)] = o

    xs_out, ys_out, zs_out, os_out = [], [], [], []
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                x = xmin + ix * step
                y = xmin + iy * step
                z = xmin + iz * step
                xs_out.append(x)
                ys_out.append(y)
                zs_out.append(z)
                os_out.append(grid[idx3(ix, iy, iz, n)])
    return xs_out, ys_out, zs_out, os_out

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
    os.makedirs(os.path.join(out_root, top + "_tra0"), exist_ok=True)
    os.makedirs(os.path.join(out_root, top + "_tra1"), exist_ok=True)

for in_path, top, sub_rel in files:
    xs, ys, zs, os_ = read_map(in_path)

    xs0, ys0, zs0, os0 = translate_and_fill(xs, ys, zs, os_, 0.0, args.xmin, args.xmax, args.step)
    out_path0 = os.path.join(out_root, top + "_tra0", sub_rel) if sub_rel else os.path.join(out_root, top + "_tra0", os.path.basename(in_path))
    write_map(out_path0, xs0, ys0, zs0, os0)

    xs1, ys1, zs1, os1 = translate_and_fill(xs, ys, zs, os_, args.step, args.xmin, args.xmax, args.step)
    out_path1 = os.path.join(out_root, top + "_tra1", sub_rel) if sub_rel else os.path.join(out_root, top + "_tra1", os.path.basename(in_path))
    write_map(out_path1, xs1, ys1, zs1, os1)

    print(f"[OK] {in_path} -> tra0(copy), tra1(+x={args.step})")

