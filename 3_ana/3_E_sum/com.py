#!/usr/bin/env python3
import sys
import os

def read_energy_file(path):
    data = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split()
            if len(t) < 4:
                continue
            rec, lig = t[0], t[1]
            try:
                E_sum_atom = float(t[-2])
                E_sum_map = float(t[-1])
            except ValueError:
                continue
            data[(rec, lig)] = (E_sum_atom, E_sum_map)
    return data

def read_limit_file(path):
    vals = {}
    with open(path, "r", encoding="utf-8") as f:
        for ln, raw in enumerate(f, 1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            t = line.split()
            if len(t) < 3:
                continue
            rec, lig = t[0], t[1]
            val = None
            for tok in reversed(t[2:]):
                try:
                    val = float(tok)
                    break
                except ValueError:
                    continue
            if val is None:
                sys.stderr.write(f"[Warn] {path}:{ln} no numeric value found\n")
                continue
            vals[(rec, lig)] = val
    return vals

def main():
    if len(sys.argv) < 3:
        print("Usage: python com.py <energy_file> <limit_file1> ... <limit_fileN> <output_file>")
        sys.exit(1)

    energy_path = sys.argv[1]
    limit_paths = sys.argv[2:-1]
    output_path = sys.argv[-1]

    if not os.path.isfile(energy_path):
        sys.stderr.write(f"[Error] Energy file not found: {energy_path}\n")
        sys.exit(1)
    for p in limit_paths:
        if not os.path.isfile(p):
            sys.stderr.write(f"[Error] Limit file not found: {p}\n")
            sys.exit(1)

    energy_data = read_energy_file(energy_path)
    limit_data_list = [read_limit_file(p) for p in limit_paths]

    # Use the parent directory name as label (e.g. max_99, max_9, ...)
    labels = [os.path.basename(os.path.dirname(p)) for p in limit_paths]

    with open(output_path, "w", encoding="utf-8") as out:
        header = (
            f"#{'receptor':20s}"
            f"{'ligand':34s}"
            f"{'E_sum(atom)':>15s}"
            f"{'E_sum(map)':>13s}"
        )
        for label in labels:
            header += f"{label:>11s}"
        out.write(header + "\n")

        for (rec, lig) in sorted(energy_data.keys()):
            E_atom, E_map = energy_data[(rec, lig)]
            line = (
                f"{rec:20s}"
                f"{lig:35s}"
                f"{E_atom:15.2f}"
                f"{E_map:13.2f}"
            )
            for vals in limit_data_list:
                v = vals.get((rec, lig), float("nan"))
                line += f"{v:11.6f}"
            out.write(line + "\n")

if __name__ == "__main__":
    main()

