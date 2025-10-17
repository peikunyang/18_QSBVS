#!/usr/bin/env python3
import sys, math, re

RECEPTOR_W = 20
LIGAND_W = 35
NUM_W = 12  # width for numeric columns

def parse_float(s: str) -> float:
    return float(s.strip())

def load_e1(path: str):
    """Return dict: (receptor, ligand_atom) -> (E_elec, E_vdw)."""
    m = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            t = line.split()
            receptor, ligand = t[0], t[1]
            e_elec, e_vdw = parse_float(t[2]), parse_float(t[3])
            m[(receptor, ligand)] = (e_elec, e_vdw)
    return m

def map_to_atom_name(name: str) -> str:
    """Convert *_tra0/tra1 (map) to *_x_0_0/_x_0_6 (atom)."""
    name = re.sub(r'_tra0$', '_x_0_0', name)
    name = re.sub(r'_tra1$', '_x_0_6', name)
    return name

def load_e2(path: str):
    """Return list of rows from E_2_map."""
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            t = line.split()
            receptor, ligand = t[0], t[1]
            e_elec, e_vdw, e_sum = parse_float(t[2]), parse_float(t[3]), parse_float(t[4])
            rows.append((receptor, ligand, e_elec, e_vdw, e_sum))
    return rows

def cap_vdw_sum(x: float) -> float:
    """Cap only VDW and SUM at 99999.99; keep NaN as-is."""
    if isinstance(x, float) and not math.isnan(x) and x > 99999.99:
        return 99999.99
    return x

def fmt_elec(x: float) -> str:
    """Format E_elec with width NUM_W and 2 decimals (no cap)."""
    return f"{x:{NUM_W}.2f}"

def fmt_vdw_sum(x: float) -> str:
    """Format VDW/SUM with cap and width NUM_W and 2 decimals."""
    return f"{cap_vdw_sum(x):{NUM_W}.2f}"

def hdr(label: str) -> str:
    """Right-align header labels to NUM_W."""
    return f"{label:>{NUM_W}s}"

def main():
    if len(sys.argv) != 4:
        print("Usage: python cmp.py <E_1_atom> <E_2_map> <output>")
        sys.exit(1)

    e1_path, e2_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    e1_map = load_e1(e1_path)
    e2_rows = load_e2(e2_path)

    with open(out_path, "w", encoding="utf-8") as w:
        # header
        w.write(
            f"#"
            f"{'receptor':<{RECEPTOR_W}s} "
            f"{'ligand':<{LIGAND_W}s} "
            f"{hdr('E_elec(atom)')} {hdr('E_elec(map)')} "
            f"{hdr('E_vdw(atom)')} {hdr('E_vdw(map)')} "
            f"{hdr('E_sum(atom)')} {hdr('E_sum(map)')}\n"
        )

        for receptor, ligand_map, e2_elec, e2_vdw, e2_sum in e2_rows:
            ligand_atom = map_to_atom_name(ligand_map)

            e1_elec, e1_vdw = (math.nan, math.nan)
            if (receptor, ligand_atom) in e1_map:
                e1_elec, e1_vdw = e1_map[(receptor, ligand_atom)]
            e1_sum = e1_elec + e1_vdw if not math.isnan(e1_elec) and not math.isnan(e1_vdw) else math.nan

            w.write(
                f"{receptor:<{RECEPTOR_W}s} "
                f"{ligand_map:<{LIGAND_W}s} "
                f"{fmt_elec(e1_elec)} {fmt_elec(e2_elec)} "
                f"{fmt_vdw_sum(e1_vdw)} {fmt_vdw_sum(e2_vdw)} "
                f"{fmt_vdw_sum(e1_sum)} {fmt_vdw_sum(e2_sum)}\n"
            )

if __name__ == "__main__":
    main()

