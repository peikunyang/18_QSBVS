import os
import argparse
import numpy as np

def read_last_column(path):
    vals = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            t = line.split()
            if len(t) >= 4:
                vals.append(float(t[-1]))
    return np.array(vals)

parser = argparse.ArgumentParser()
parser.add_argument("--new_root", required=True)
parser.add_argument("--old_root", required=True)
parser.add_argument("--tol", type=float, default=1e-8)
args = parser.parse_args()

new_root = os.path.abspath(args.new_root)
old_root = os.path.abspath(args.old_root)

mapping = {
    "_tra0": "_x_0_0",
    "_tra1": "_x_0_6"
}

dirs = sorted([d for d in os.listdir(new_root) if os.path.isdir(os.path.join(new_root, d))])

for d in dirs:
    match = None
    for k, v in mapping.items():
        if k in d:
            match = d.replace(k, v)
            break
    if match is None:
        print(f"[SKIP] {d} has no matching pattern")
        continue

    new_dir = os.path.join(new_root, d)
    old_dir = os.path.join(old_root, match)
    if not os.path.exists(old_dir):
        print(f"[MISS] no old dir for {d}")
        continue

    files = sorted(os.listdir(new_dir))
    for f in files:
        new_path = os.path.join(new_dir, f)
        old_path = os.path.join(old_dir, f)
        if not os.path.exists(old_path):
            print(f"[MISS] {f} missing in old {match}")
            continue

        a = read_last_column(new_path)
        b = read_last_column(old_path)
        if len(a) != len(b):
            print(f"[DIFF] {d}/{f} length {len(a)} vs {len(b)}")
            continue
        diff = np.abs(a - b)
        maxdiff = diff.max() if len(diff) > 0 else 0.0
        if maxdiff > args.tol:
            idx = int(np.argmax(diff))
            print(f"[DIFF] {d}/{f} maxdiff={maxdiff:.3e} at line {idx+1}")
        else:
            print(f"[OK]   {d}/{f}")

