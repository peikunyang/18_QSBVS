import os

shifts = [0.0, 0.6]

in_dir = "../2_rot/pdbqt"
out_dir = "pdbqt"
os.makedirs(out_dir, exist_ok=True)

for fn in os.listdir(in_dir):
    in_path = os.path.join(in_dir, fn)
    if not os.path.isfile(in_path):
        continue
    if fn.endswith("_x0") or fn.endswith("_x0p6"):
        continue

    with open(in_path, "r") as fr:
        lines = fr.readlines()

    for dx in shifts:
        suffix = "_0" if dx == 0.0 else "_1"
        out_path = os.path.join(out_dir, fn + suffix)

        with open(out_path, "w") as fw:
            for line in lines:
                if not line.strip():
                    fw.write(line)
                    continue

                p = line.split()
                if len(p) < 7:
                    fw.write(line)
                    continue

                typ = p[0]
                x = float(p[1]) + dx
                y = float(p[2])
                z = float(p[3])
                q = float(p[4])
                r_str = p[5]
                e_str = p[6]

                fw.write(
                    f"{typ:>3} "
                    f"{x:12.3f} {y:12.3f} {z:12.3f} "
                    f"{q:9.3f} {r_str:>7} {e_str:>12}\n"
                )

