#!/usr/bin/env python3
import numpy as np
import argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", help="輸入檔案路徑")
    args = ap.parse_args()

    # 讀檔
    data = []
    with open(args.infile, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                data.append(float(line))
            except ValueError:
                pass  # 如果遇到非數值的行就跳過

    arr = np.array(data, dtype=float)

    # 計算 norm
    norm = np.linalg.norm(arr)

    print(f"Norm = {norm}")

if __name__ == "__main__":
    main()

