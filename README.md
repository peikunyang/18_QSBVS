This repository contains the data preparation, energy evaluation, and analysis pipeline used in the study  
**“Map-Based Quantum Inner-Product Scoring for Structure-Based Virtual Screening.”**

The workflow implements a map-based reformulation of protein–ligand interaction energies, including electrostatic and van der Waals contributions, and supports both classical and quantum inner-product–based energy evaluation.

---

## Repository Structure

The repository is organized into three main stages:

```
18_QSBVS/
├── 1_data/      # Input structures, force-field data, and grid maps
├── 2_energy/    # Electrostatic and van der Waals energy calculations
├── 3_ana/       # Data extraction and analysis
```

Each stage corresponds to a distinct step in the structure-based virtual screening pipeline.

---

## 1. Data Preparation (`1_data/`)

This directory contains all structural inputs and precomputed grid data required for energy evaluation.

### `1_data/3_pdbqt/4_cen/`
- Contains centered **PDBQT** files.
- Includes:
  - **2 protein structures**
  - **4 ligand structures**
- All structures are aligned and centered in a common Cartesian reference frame.

### `1_data/5_map/3_pdbqt/`
- Contains files required for map-based energy evaluation:
  - Protein and ligand coordinate files
  - Force-field parameter files
  - Precomputed grid maps
- Includes data for:
  - **2 proteins**
  - **4 ligands**
- These maps are used to evaluate electrostatic and van der Waals interactions on a shared grid.

---

## 2. Energy Evaluation (`2_energy/`)

This directory contains scripts for computing interaction energies using the map-based formulation.

- Electrostatic interaction energy
- van der Waals interaction energy
- Energies are evaluated as inner products between:
  - Receptor potential maps
  - Ligand charge and atom-type occupancy vectors

Both probability-based (ideal) and finite-shot sampling–based evaluations are supported.

---

## 3. Data Extraction and Analysis (`3_ana/`)

This directory contains scripts for post-processing and analysis of computed energies, including:

- Extraction of energy values from raw outputs
- Statistical analysis under different truncation thresholds
- Evaluation of finite-sampling effects
- Ranking and comparison of receptor–ligand configurations

## Notes

- The repository is designed to mirror the computational pipeline described in the associated manuscript.
- Quantum inner-product estimation is implemented using the Hadamard test, as described in the Supplementary Information of the paper.
- This repository serves as a reference implementation and data archive rather than a standalone docking package.

---

## Data and Software Availability

All data and scripts supporting the study are available in this repository:

👉 https://github.com/peikunyang/18_QSBVS

---

## Author

**Pei-Kun Yang**  
Email: peikun@isu.edu.tw  
Alternate: peikun6416@gmail.com
