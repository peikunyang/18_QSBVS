#!/usr/bin/env python3
import os
import sys
from collections import defaultdict

def read_third_from_last_column(file_path):
    vals = []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            # Skip empty, comment, or separator lines
            if not s or s.startswith("#") or set(s) <= {"-", " "}:
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            try:
                vals.append(float(parts[-3]))  # third from last column
            except ValueError:
                continue
    return vals

def parse_name(fname):
    base = os.path.splitext(fname)[0]
    parts = base.split("_")
    if len(parts) == 2 and parts[0] == "energy":
        try:
            g = int(parts[1])
            return (g, None)
        except:
            return (999999, None)
    if len(parts) == 3 and parts[0] == "energy":
        try:
            g = int(parts[1])
            i = int(parts[2])
            return (g, i)
        except:
            return (999999, None)
    return (999999, None)

def average_rank(vals, target):
    sorted_vals = sorted(vals)
    idxs = [i for i, v in enumerate(sorted_vals) if v == target]
    if not idxs:
        return None
    return (idxs[0] + idxs[-1]) / 2 + 1  # 1-based average rank

def main():
    if len(sys.argv) < 3:
        print("Usage: python rank_energy.py <input_dir> <output_file>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2]

    files = []
    for f in os.listdir(input_dir):
        if not f.startswith("energy"):
            continue
        path = os.path.join(input_dir, f)
        if os.path.isfile(path):
            files.append(f)
    files.sort(key=lambda x: parse_name(x))

    group_map = defaultdict(list)
    for f in files:
        g, i = parse_name(f)
        if g == 999999:
            continue
        group_map[g].append((i, f))

    with open(output_file, "w", encoding="utf-8") as out:
        if 0 in group_map:
            for i, fname in sorted(group_map[0], key=lambda t: (t[0] is None, t[0] if t[0] is not None else -1)):
                path = os.path.join(input_dir, fname)
                vals = read_third_from_last_column(path)
                if not vals:
                    continue
                target = vals[0]
                rank = average_rank(vals, target)
                if rank is None:
                    continue
                out.write(f"energy_0 {rank:5.2f}\n")
            del group_map[0]

        for g in sorted(group_map.keys()):
            ranks = []
            ordered = sorted(group_map[g], key=lambda t: (t[0] is None, t[0] if t[0] is not None else -1))
            for i, fname in ordered:
                path = os.path.join(input_dir, fname)
                vals = read_third_from_last_column(path)
                if not vals:
                    continue
                target = vals[0]
                rank = average_rank(vals, target)
                if rank is None:
                    continue
                ranks.append((i, rank))
            if not ranks:
                continue
            ranks_sorted = sorted([r for r in ranks if r[0] is not None], key=lambda t: t[0])
            rank_vals = [r[1] for r in ranks_sorted]
            if not rank_vals:
                continue
            avg = sum(rank_vals) / len(rank_vals)
            parts = " ".join(f"{v:5.2f}" for v in rank_vals)
            out.write(f"energy_{g}  {avg:6.2f} {parts}\n")

if __name__ == "__main__":
    main()

