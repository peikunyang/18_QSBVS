#!/usr/bin/env python3
import os
import numpy as np

base_dir = "../1_data"
X_list = ["max_99", "max_9", "max_1", "max_0.4"]
Y_list = ["R_5_matrix_mul", "R_3_Hadamard_test", "R_4_SWAP_test"]
output_file = "energy"

def extract_energies(path):
    data = {}
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            label = parts[0]
            if label.startswith("energy_"):
                idx = int(label.split("_")[1])
                if 1 <= idx <= 7:
                    vals = [float(x) for x in parts[1:]]
                    data[label] = np.mean(vals)
    return [data.get(f"energy_{i}", np.nan) for i in range(1, 8)]

with open(output_file, "w") as out:
    out.write(f"{'Method':<22}{'Threshold':>10}{'energy_1':>12}{'energy_2':>12}{'energy_3':>12}{'energy_4':>12}{'energy_5':>12}{'energy_6':>12}{'energy_7':>12}\n")
    for Y in Y_list:
        for X in X_list:
            file_path = os.path.join(base_dir, X, Y)
            if not os.path.exists(file_path):
                continue
            energies = extract_energies(file_path)
            out.write(f"{Y:<22}{X.replace('max_',''):>10}" + "".join([f"{v:12.2f}" for v in energies]) + "\n")

