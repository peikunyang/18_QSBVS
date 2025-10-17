# 🧬 Map-Based Quantum Inner-Product Scoring for Structure-Based Virtual Screening

This repository contains the data structure, preprocessing workflow, and analysis scripts used in the study:

> **Yang, P.-K. (2025).** *Map-Based Quantum Inner-Product Scoring for Structure-Based Virtual Screening.*

---

## 📂 Repository Structure

```
.
├── 1_data/
│   ├── 1_ori/
│   ├── 2_MD/
│   └── 3_pdbqt/
├── 2_energy/
│   ├── 1_energy/
│   └── 2_limit_energy/
└── 3_ana/
```

---

### 🧩 `1_data/`

#### **1_ori/**
- Contains the original **protein structure (PDB ID: 8F76)** used as the receptor template.

#### **2_MD/**
- Includes **molecular dynamics (MD) simulations** performed using CHARMM-GUI and OpenMM.
- Two systems were generated:
  - Protein–ligand complex.
  - Apo protein without ligand.
- Files are **not included** in the repository due to storage constraints.  
  → Contact the author to obtain them.

#### **3_pdbqt/**
- Contains `.pdbqt` files and corresponding map files generated for AutoDock-based energy evaluation:
  - **2 protein conformations** (`.pdbqt`)
  - **2 ligands**, each with **2 conformations** → total of **4 ligand `.pdbqt` files**
  - **Receptor maps** (derived from 2 protein conformations)
  - **Ligand maps** (4 ligand conformations)
  - After applying **rotations** and **translations**, a total of **16 ligand map configurations** were generated.

---

### ⚡ `2_energy/`

This directory contains classical and quantum-derived energy evaluations.

#### **1_energy/**
- **1_atom/** — Energy calculated using **atom-based pairwise summation**.
- **2_map/** — Energy calculated using **grid map-based dot products**.

#### **2_limit_energy/**
Each subdirectory corresponds to a **truncation threshold** for potential normalization (`|Umax|`, `|Nmax|`):

```
max_99/
max_9/
max_1/
max_0.4/
```

Each folder includes energy results computed by multiple methods:

| Folder | Description |
|--------|--------------|
| `1_limit_max` | Classical map-based energy (reference) |
| `2_dot_product` | Inner-product (dot product) method |
| `3_Hadamard_test` | Quantum estimation using **Hadamard test** |
| `4_SWAP_test` | Quantum estimation using **SWAP test** |
| `5_matrix_mul` | Quantum estimation using **GEMM (matrix multiplication)** |

These results correspond to the data shown in Tables II and III of the paper, including convergence and ranking analyses across thresholds.

---

### 📊 `3_ana/`

Contains post-processing and **analysis scripts** for:

- Statistical comparison of energy evaluation methods  
- Ranking-based performance metrics (mean rank, shot efficiency)  
- Visualization and figure generation for the manuscript  

---

## ⚙️ Computational Workflow Overview

1. **Protein & Ligand Preparation**  
   - CHARMM-GUI and OpenMM simulations using CHARMM36m and CGenFF.  
   - Two receptor conformations and two ligands × two conformations each.

2. **Grid Map Construction**  
   - AutoDock4 maps (electrostatic and vdW) on an 8×8×8 grid (512 points).  
   - Ligands rotated and translated to generate 16 spatial configurations.

3. **Energy Evaluation Methods**
   - **Classical:** Atom-based and map-based dot product.  
   - **Quantum:**  
     - **GEMM Method**  
     - **Hadamard Test**  
     - **SWAP Test**

4. **Shot-based Convergence Analysis**  
   - Measurements from 10¹ to 10⁷ shots.  
   - Mean rank used to evaluate identification accuracy of lowest-energy configurations.

---

## 📘 Reference Paper

> **Yang, P.-K. (2025).** *Map-Based Quantum Inner-Product Scoring for Structure-Based Virtual Screening.*  
>  
> Keywords: AutoDock grid maps, Hadamard test, SWAP test, GEMM  
>  
> **Contact:** peikun@isu.edu.tw

---

## 🧠 Citation

If you use this repository or its methods, please cite:

```
@article{Yang2025_QSBVS,
  author = {Pei-Kun Yang},
  title = {Map-Based Quantum Inner-Product Scoring for Structure-Based Virtual Screening},
  year = {2025},
  keywords = {AutoDock, Hadamard Test, SWAP Test, GEMM, Quantum Computing, Virtual Screening}
}
```

---

## 📬 Contact

For questions or access to full MD trajectory files:  
📧 **peikun@isu.edu.tw**
