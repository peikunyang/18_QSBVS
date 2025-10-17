import os, sys, torch, numpy as np

N_qubit = 11
DIM = 2**N_qubit

def Ry_gate(theta):
    c = torch.cos(theta/2); s = torch.sin(theta/2)
    return torch.tensor([[c, -s], [s, c]], dtype=torch.float64)

def I_gate():
    return torch.tensor([[1., 0.], [0., 1.]], dtype=torch.float64)

def recursive_amplitude_encoding(states, alpha, code):
    length = len(states)
    if length == 2:
        alpha.append(-2.0 * np.arctan2(states[1], states[0]))
        code.append(int(N_qubit - np.log2(length)))
    elif length > 2:
        half = length // 2
        norm_top = np.linalg.norm(states[:half])
        norm_bottom = np.linalg.norm(states[half:])
        alpha.append(-2.0 * np.arctan2(norm_bottom, norm_top))
        code.append(int(N_qubit - np.log2(length)))
        recursive_amplitude_encoding(states[:half], alpha, code)
        recursive_amplitude_encoding(states[half:], alpha, code)

def gen_angles(states):
    alpha, code = [], []
    recursive_amplitude_encoding(states, alpha, code)
    return alpha, code

def param_con(p, code):
    params = []
    for j in range(N_qubit):
        par = []
        for i in range(len(code)):
            if code[i] == j:
                par.append(p[i])
        params.append(torch.tensor(par, dtype=torch.float64))
    return params

def read_states(path):
    v = np.loadtxt(path, dtype=np.float64)
    if v.size != DIM: raise ValueError(f"len {v.size} != {DIM}: {path}")
    nrm = np.linalg.norm(v)
    if nrm == 0: raise ValueError(f"zero vector: {path}")
    return v / nrm, nrm

def Tensor_Prod(ry, lay):
    out = ry
    for _ in range(lay):
        out = torch.kron(out, I_gate())
    return out

def Circuit(par):
    ry = Ry_gate(par[0][0])
    qc1 = Tensor_Prod(ry, N_qubit-1)
    for i in range(1, N_qubit):
        qc2 = torch.zeros((DIM, DIM), dtype=torch.float64)
        m = par[i].shape[0]
        lay = int(N_qubit - np.log2(m))
        for j in range(m):
            blk = Tensor_Prod(Ry_gate(par[i][j]), lay-1)
            a = j * (2**lay)
            b = (j + 1) * (2**lay)
            qc2[a:b, a:b] = blk
        qc1 = qc1 @ qc2
    return qc1

def process_all(input_dir, output_dir, tol=1e-6):
    os.makedirs(output_dir, exist_ok=True)
    files = sorted([f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))])
    for fname in files:
        fpath = os.path.join(input_dir, fname)
        u, norm = read_states(fpath)
        angles, code = gen_angles(u)
        pars = param_con(angles, code)
        U = Circuit(pars)
        diff_row = np.linalg.norm(U[0, :].numpy() - u)
        out_name = os.path.splitext(fname)[0] + ".pt"
        torch.save(U, os.path.join(output_dir, out_name))
        print(f"{fname}: input norm={norm:.6f}, first-row diff={diff_row:.3e}")
    print("All matrices generated successfully.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python unitary_origrow.py <input_dir> <output_dir>"); sys.exit(1)
    process_all(sys.argv[1], sys.argv[2])

