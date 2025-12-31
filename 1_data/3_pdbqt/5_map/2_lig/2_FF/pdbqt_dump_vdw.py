#!/usr/bin/env python3
import os, glob, argparse, sys
from collections import defaultdict, Counter

def load_params(path):
    params = {}
    if path and os.path.isfile(path):
        with open(path, "r", encoding="latin-1") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith(("#", "!", ";")):
                    continue
                t = s.split()
                if t[0].lower().startswith("atom_par") and len(t) >= 4:
                    try:
                        at = t[1]
                        rii = float(t[2])
                        eps = float(t[3])
                        params[at] = (rii, eps)
                    except:
                        pass
    if not params:
        params.update({"C": (4.00, 0.150), "OA": (3.20, 0.200), "HD": (2.00, 0.020)})
    return params

def parse_atom_line(line):
    rec = line[:6].strip()
    if rec not in ("ATOM", "HETATM"):
        return None
    x = y = z = None
    try:
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
    except:
        t = line.split()
        if len(t) >= 12:
            x = float(t[6]); y = float(t[7]); z = float(t[8])
        else:
            return None
    t = line.split()
    typ = t[-1]
    q = None
    try:
        q = float(t[-2])
    except:
        q = float("nan")
    return typ, x, y, z, q

def process_file(inpath, outpath, params):
    rows = []
    missing_counter = Counter()
    with open(inpath, "r", encoding="latin-1") as f:
        for line in f:
            p = parse_atom_line(line)
            if not p:
                continue
            typ, x, y, z, q = p
            R, eps = params.get(typ, (None, None))
            if R is None or eps is None:
                missing_counter[typ] += 1
            rows.append((typ, x, y, z, q, R, eps))
    with open(outpath, "w", encoding="utf-8") as w:
        for typ, x, y, z, q, R, eps in rows:
            r_str = "" if R is None else f"{R:7.2f}"
            e_str = "" if eps is None else f"{eps:12.3f}"
            w.write(f"{typ:>3} {x:12.3f} {y:12.3f} {z:12.3f} {q:9.3f} {r_str:>7} {e_str:>12}\n")
    return missing_counter

ap = argparse.ArgumentParser()
ap.add_argument("--indir", required=True)
ap.add_argument("--outdir", default="vdw")
ap.add_argument("--param", default="")
ap.add_argument("--summary", default="missing_vdw_types.txt")
args = ap.parse_args()

os.makedirs(args.outdir, exist_ok=True)
params = load_params(args.param)
files = sorted(glob.glob(os.path.join(args.indir, "*.pdbqt")))

global_missing = defaultdict(lambda: {"count": 0, "files": set()})

for fp in files:
    bn = os.path.basename(fp)
    outp = os.path.join(args.outdir, os.path.splitext(bn)[0] + "_FF")
    mcnt = process_file(fp, outp, params)
    if mcnt:
        sys.stderr.write("[WARN] " + bn + ": missing vdW params for " +
                         ", ".join(f"{t}({c})" for t, c in mcnt.items()) + "\n")
        for t, c in mcnt.items():
            global_missing[t]["count"] += c
            global_missing[t]["files"].add(bn)

if global_missing:
    with open(args.summary, "w", encoding="utf-8") as s:
        s.write("# type\ttotal_missing\tfiles\n")
        for t, info in sorted(global_missing.items()):
            s.write(f"{t}\t{info['count']}\t{','.join(sorted(info['files']))}\n")
    sys.stderr.write("[WARN] Summary -> " + args.summary + "\n")

print(f"Done: {len(files)} files -> {args.outdir}")

