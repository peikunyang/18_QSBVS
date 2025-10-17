#!/usr/bin/env python3
import os
import sys
import re
from math import isnan

def read_second_column(path):
    vals = []
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                if not line.strip():
                    continue
                t = line.strip().split()
                if len(t) < 2:
                    continue
                try:
                    v = float(t[1])
                except ValueError:
                    continue
                vals.append(v)
    except FileNotFoundError:
        pass
    return vals

def discover_max_dirs(base_dir):
    """Return a sorted list of max_* directory names (e.g., ['max_99','max_9','max_1','max_0.4'])."""
    dirs = []
    if not os.path.isdir(base_dir):
        return dirs
    for name in os.listdir(base_dir):
        full = os.path.join(base_dir, name)
        if not os.path.isdir(full):
            continue
        if not name.startswith("max_"):
            continue
        # extract numeric part after "max_"
        num_str = name[4:]
        # allow forms like "0.4" or "99"
        try:
            num_val = float(num_str)
        except ValueError:
            continue
        dirs.append((name, num_val))
    # sort by numeric value descending (e.g., 99 > 9 > 1 > 0.4)
    dirs.sort(key=lambda x: x[1], reverse=True)
    return [name for name, _ in dirs]

def main():
    if len(sys.argv) < 3:
        print("Usage: python comb.py <subdir> <output_file>", file=sys.stderr)
        sys.exit(1)

    subdir = sys.argv[1]
    output_file = sys.argv[2]

    base_dir = "../../../2_energy/2_limit_energy"

    # auto-discover existing max_* directories
    max_dirs = discover_max_dirs(base_dir)
    if not max_dirs:
        print(f"[Error] No max_* directories found under: {base_dir}", file=sys.stderr)
        sys.exit(1)

    all_vals = []
    for mdir in max_dirs:
        rank_path = os.path.join(base_dir, mdir, "2_energy", subdir, "ana", "rank")
        if not os.path.exists(rank_path):
            print(f"[Skip] {rank_path}", file=sys.stderr)
            all_vals.append([])
            continue
        vals = read_second_column(rank_path)
        all_vals.append(vals)

    num_rows = max((len(v) for v in all_vals if v), default=0)

    # formatting
    lead_w = 6
    col_w = 10

    with open(output_file, "w", encoding="utf-8") as out:
        # header
        out.write(f"{'':<{lead_w}s}" + "".join(f"{mdir:>{col_w}s}" for mdir in max_dirs) + "\n")
        # first data row label "pro"
        out.write(f"{'pro':<{lead_w}s}")
        for vals in all_vals:
            if len(vals) > 0:
                out.write(f"{vals[0]:{col_w}.2f}")
            else:
                out.write(f"{'nan':>{col_w}s}")
        out.write("\n")

        # subsequent rows indexed from 1
        for i in range(1, num_rows):
            out.write(f"{i:<{lead_w}d}")
            for vals in all_vals:
                if i < len(vals):
                    out.write(f"{vals[i]:{col_w}.2f}")
                else:
                    out.write(f"{'nan':>{col_w}s}")
            out.write("\n")

    print(f"Done -> {output_file}")

if __name__ == "__main__":
    main()

