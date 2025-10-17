#!/usr/bin/env python3
import os
import sys
import math

MAP_FILES = ["map_elec", "map_vdw_C", "map_vdw_OA"]

def read_last_column_512(path, limit):
    vals = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if len(vals) >= 512:
                break
            t = line.strip().split()
            if not t:
                continue
            try:
                v = float(t[-1])
            except ValueError:
                continue
            # 截斷
            if v > limit:
                v = limit
            if v < -limit:
                v = -limit
            vals.append(v)
    if len(vals) < 512:
        vals.extend([0.0] * (512 - len(vals)))
    return vals[:512]

def main():
    if len(sys.argv) < 2:
        print("Usage: python comb.py <input_dir> [limit]")
        sys.exit(1)

    input_dir = sys.argv[1]
    limit = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4  # ← 預設0.4
    print(f"Using cutoff limit: ±{limit}")

    output_dir = "./maps"
    os.makedirs(output_dir, exist_ok=True)

    summary_lines = []
    for sub in sorted(os.listdir(input_dir)):
        sub_path = os.path.join(input_dir, sub)
        if not os.path.isdir(sub_path):
            continue

        values = []
        ok = True
        for fname in MAP_FILES:
            fpath = os.path.join(sub_path, fname)
            if not os.path.exists(fpath):
                ok = False
                break
            values.extend(read_last_column_512(fpath, limit))
        if not ok:
            continue

        # add 512 zeros
        values.extend([0.0] * 512)

        # write combined file (2048 lines)
        out_path = os.path.join(output_dir, sub)
        with open(out_path, "w", encoding="utf-8") as f:
            for v in values:
                f.write(f"{v:10.6f}\n")

        # sqrt of sum of squares
        ss = sum(v * v for v in values)
        rss = math.sqrt(ss)
        summary_lines.append(f"{sub:<24}{rss:16.6f}")

    with open("lengths", "w", encoding="utf-8") as f:
        for line in summary_lines:
            f.write(line + "\n")

if __name__ == "__main__":
    main()

