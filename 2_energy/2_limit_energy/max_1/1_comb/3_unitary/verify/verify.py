import os
import sys
import torch
import numpy as np

N = 2048

def is_unitary(U, tol=1e-6):
    if U.shape != (N, N):
        return False
    I = torch.eye(N, dtype=U.dtype)
    err = torch.norm(torch.mm(U.T, U) - I)
    return err.item() < tol

def check_all(unitary_dir, vector_dir):
    files = sorted([f for f in os.listdir(unitary_dir) if f.endswith(".pt")])
    for f in files:
        base = os.path.splitext(f)[0]
        u_path = os.path.join(unitary_dir, f)
        v_path = os.path.join(vector_dir, base)
        if not os.path.exists(v_path):
            print(f"{f}: missing original file ({v_path})")
            continue
        U = torch.load(u_path, map_location="cpu", weights_only=True)
        if not isinstance(U, torch.Tensor):
            print(f"{f}: not tensor")
            continue
        ok = is_unitary(U)
        v = np.loadtxt(v_path, dtype=np.float64)  # ← 不再正規化
        first_row = U[0, :].numpy()
        diff = np.linalg.norm(first_row - v)
        same = diff < 1e-6
        if ok and same:
            print(f"{f}: OK")
        elif not ok:
            print(f"{f}: NOT unitary ({tuple(U.shape)})")
        else:
            print(f"{f}: first row mismatch (diff={diff:.3e})")
    print("Check complete.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python verify.py <unitary_dir> <vector_dir>")
        sys.exit(1)
    check_all(sys.argv[1], sys.argv[2])

