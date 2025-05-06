# Protein Surface Charge Calculator

This Python script (`Protein_surface_charge.py`) ranks protein models by their absolute net surface charge. It accepts either mmCIF (`.cif`) or PDB (`.pdb`) files as input, uses APBS when available, and falls back to a simple residue-based estimate if APBS fails.

---

## Features

- **Supports `.cif` and `.pdb`** input files
- **CIF → PDB conversion** via Biopython
- **PDB → PQR conversion** using `pdb2pqr`
- **Electrostatic calculation** with APBS
- **Residue-based fallback** when APBS or conversion fails
- **Automated ranking** of multiple models by |net surface charge|

---

## Requirements

- **Python 3.7+**
- **Biopython**
- **pdb2pqr** (with PARSE force field)
- **APBS** (version ≥ 3.0, optional but recommended)
- Unix-like shell (`bash`/`zsh`) for environment variables

---

## Installation

1. **Clone or copy** this script into your working directory.

2. **Create & activate** a conda environment:

   ```bash
   conda create -n surfcharge python=3.9 -y
   conda activate surfcharge
   ```

3. **Install Python package**:

   ```bash
   pip install biopython
   ```

4. **Install `pdb2pqr`**:

   ```bash
   conda install -c conda-forge pdb2pqr -y
   ```

5. **Install APBS** (optional):

   ```bash
   # macOS Homebrew
   brew install apbs

   # or conda-forge
   conda install -c conda-forge apbs -y
   ```

6. **(Homebrew only)** add APBS tools to your PATH:

   ```bash
   export APBS_TOOLS_HOME=/opt/homebrew/opt/apbs/share/apbs/tools
   export PATH="$PATH:$APBS_TOOLS_HOME/bin"
   ```

7. **Verify**:

   ```bash
   pdb2pqr --version
   apbs --version   # if installed
   python - <<<'import Bio; print(Bio.__version__)'
   ```

---

## Usage

1. **Place** your `.cif` and/or `.pdb` files in the same folder as `Protein_surface_charge.py`.

2. **Run**:

   ```bash
   python Protein_surface_charge.py
   ```

3. **Output**:

   - A `converted_pdb/` directory containing intermediate `.pdb` files from any `.cif` inputs.
   - `.pqr` and temporary APBS input/output files in the working directory.
   - Printed ranking by absolute net surface charge, with method tag (`APBS` or `residue-count`).

---

## Interpreting Results

Each line in the output shows:

```
1. filename.ext    Charge: +12.34  (APBS)
```

- **Charge**: net surface charge in elementary charge units (e);
- **Method**: the calculation source:
  - `APBS`: high‑accuracy Poisson–Boltzmann solution
  - `residue-count`: quick estimate from charged residues at pH ≈ 7

Models are sorted by descending |net charge|.

---

## Drawbacks & Caveats

- **Hard‑coded APBS settings** (grid size, boundaries) may need tuning for large/complex systems.
- **Installation complexity**: requires external binaries (`pdb2pqr`, optionally `apbs`).
- **Residue‑count fallback** ignores pKa shifts, local environments, and partial protonation.
- **Performance**: APBS runs can be slow, especially for many or large proteins.

---

## License

MIT © Your Name
