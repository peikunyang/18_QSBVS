#!/usr/bin/env python3
import os
import sys

def read_second_column(path):
    vals = []
    with open(path, "r") as f:
        for line in f:
            t = line.strip().split()
            if len(t) < 2:
                continue
            try:
                v = float(t[1])
            except ValueError:
                continue
            vals.append(v)
    return vals

def main():
    if len(sys.argv) < 3:
        print("Usage: python comb.py <subdir> <output_file>")
        sys.exit(1)

    subdir = sys.argv[1]
    output_file = sys.argv[2]

    base_dir = "../../../2_energy/2_limit_energy"
    max_dirs = ["max_99", "max_9", "max_1", "max_0.4", "max_0.1"]

    all_vals = []
    for mdir in max_dirs:
        rank_path = os.path.join(base_dir, mdir, "2_energy", subdir, "ana", "rank")
        if not os.path.exists(rank_path):
            print(f"[Skip] {rank_path}")
            all_vals.append([])
            continue
        vals = read_second_column(rank_path)
        all_vals.append(vals)

    num_rows = max((len(v) for v in all_vals if v), default=0)

    with open(output_file, "w") as out:
        # 第一列：空白 + max 標題
        out.write(f"{'':<6s}" + "".join(f"{mdir:>10s}" for mdir in max_dirs) + "\n")
        # 第二列：pro
        out.write(f"{'pro':<6s}")
        for vals in all_vals:
            if len(vals) > 0:
                out.write(f"{vals[0]:10.2f}")
            else:
                out.write(f"{'nan':>10s}")
        out.write("\n")
        # 從第2個值開始逐列輸出
        for i in range(1, num_rows):
            out.write(f"{i:<6d}")
            for vals in all_vals:
                if i < len(vals):
                    out.write(f"{vals[i]:10.2f}")
                else:
                    out.write(f"{'nan':>10s}")
            out.write("\n")

    print(f"Done -> {output_file}")

if __name__ == "__main__":
    main()

