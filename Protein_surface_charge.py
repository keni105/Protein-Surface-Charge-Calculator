#!/usr/bin/env python3
"""
Rank protein models by absolute net surface charge.
Uses APBS when available, falling back to residue-based count if APBS fails.
Supports either .cif or .pdb inputs.
"""
import os
import subprocess
import glob
from collections import Counter
from Bio.PDB import MMCIFParser, PDBIO, PDBParser

def convert_cif_to_pdb(cif_file, output_dir="converted_pdb"):
    """
    Convert an mmCIF file to PDB format using Biopython.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        parser = MMCIFParser(QUIET=True)
        struct_id = os.path.splitext(os.path.basename(cif_file))[0]
        structure = parser.get_structure(struct_id, cif_file)
        pdb_path = os.path.join(output_dir, f"{struct_id}.pdb")
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)
        return pdb_path
    except Exception as e:
        print(f"✗  CIF→PDB failed for {cif_file}: {e}")
        return None

def fallback_charge(pdb_file):
    """
    Simple residue-based net charge estimate at pH≈7.
    Arg/Lys = +1, Asp/Glu = -1, His ≈ +0.1
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_file)
    counts = Counter()
    for res in structure.get_residues():
        # only standard residues
        if res.id[0] == ' ':
            counts[res.resname] += 1
    charge = (
        counts.get('ARG', 0) * 1.0 +
        counts.get('LYS', 0) * 1.0 +
        counts.get('ASP', 0) * -1.0 +
        counts.get('GLU', 0) * -1.0 +
        counts.get('HIS', 0) * 0.1
    )
    return charge

def calculate_surface_charge(pdb_file):
    """
    Calculate surface charge: first try APBS, else fallback.
    Returns (charge, method_str).
    """
    base = os.path.splitext(os.path.basename(pdb_file))[0]
    pqr = f"{base}.pqr"
    try:
        # PDB→PQR
        subprocess.run(
            ["pdb2pqr", "--ff=PARSE", pdb_file, pqr],
            check=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        # Prepare minimal APBS input
        in_file = f"{base}.in"
        with open(in_file, 'w') as inp:
            inp.write(f"""
read
    mol pqr {pqr}
end
energy
    calculate total charge
end
quit
""")
        # Run APBS
        res = subprocess.run(
            ["apbs", in_file],
            capture_output=True, text=True
        )
        # Parse net‐charge line
        for line in res.stdout.splitlines():
            if "Net charge" in line:
                # e.g. '  Net charge -1.30e+01 e'
                parts = line.split()
                # second‐to‐last token is the number
                charge = float(parts[-2])
                return charge, 'APBS'
        raise RuntimeError("APBS did not report net charge")
    except Exception:
        # fallback
        print(f"→ Using fallback residue-based estimate for {pdb_file}")
        return fallback_charge(pdb_file), 'residue-count'

def rank_proteins_by_charge(files):
    """
    Accepts a list of .cif or .pdb files.
    Returns list of (orig_filename, charge, method) sorted by |charge| desc.
    """
    results = []
    for fn in files:
        print(f"→ Processing {fn}")
        ext = os.path.splitext(fn)[1].lower()
        if ext == ".cif":
            pdb = convert_cif_to_pdb(fn)
            if not pdb:
                continue
        elif ext == ".pdb":
            pdb = fn
        else:
            print(f"✗  Skipping unsupported file {fn}")
            continue
        charge, method = calculate_surface_charge(pdb)
        results.append((os.path.basename(fn), charge, method))
    # sort by absolute charge descending
    results.sort(key=lambda x: abs(x[1]), reverse=True)
    return results

if __name__ == '__main__':
    # gather both .cif and .pdb
    files = sorted(glob.glob("*.cif") + glob.glob("*.pdb"))
    if not files:
        print("No .cif or .pdb files found in current directory.")
        exit(1)

    ranked = rank_proteins_by_charge(files)
    print("\nRanked Proteins by Surface Charge:")
    for i, (name, charge, method) in enumerate(ranked, 1):
        print(f"{i}. {name:<30} Charge: {charge:>7.2f}  ({method})")

    # Write to text file
    out_txt = "surface_charge_ranking.txt"
    with open(out_txt, "w") as outf:
        outf.write("Rank\tModel\tCharge\tMethod\n")
        for i, (name, q, method) in enumerate(ranked, 1):
            outf.write(f"{i}\t{name}\t{q:.2f}\t{method}\n")

    print(f"\nResults written to {out_txt}")
