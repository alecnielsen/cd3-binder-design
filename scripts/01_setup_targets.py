#!/usr/bin/env python3
"""Step 1: Setup target structures.

Downloads and prepares CD3 target structures for design:
- 1XIW: CD3εδ heterodimer
- 1SY6: CD3εγ + OKT3 Fab complex

Usage:
    python scripts/01_setup_targets.py [--output-dir data/targets]
"""

import argparse
from pathlib import Path
import urllib.request


# PDB download URL template
PDB_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.pdb"


def download_pdb(pdb_id: str, output_path: Path) -> bool:
    """Download PDB file from RCSB.

    Args:
        pdb_id: PDB identifier (e.g., "1XIW").
        output_path: Path to save the file.

    Returns:
        True if successful.
    """
    url = PDB_URL_TEMPLATE.format(pdb_id=pdb_id)

    try:
        print(f"  Downloading {pdb_id} from {url}...")
        urllib.request.urlretrieve(url, output_path)
        print(f"  Saved to: {output_path}")
        return True
    except Exception as e:
        print(f"  ERROR: Failed to download {pdb_id}: {e}")
        return False


def extract_cd3_epsilon(pdb_path: Path, output_path: Path, chain_id: str = "E") -> bool:
    """Extract CD3ε chain from PDB file.

    Args:
        pdb_path: Path to input PDB.
        output_path: Path to save extracted chain.
        chain_id: Chain ID for CD3ε.

    Returns:
        True if successful.
    """
    try:
        with open(pdb_path, "r") as f:
            lines = f.readlines()

        # Filter to only include the specified chain
        extracted_lines = []
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[21] == chain_id:
                    extracted_lines.append(line)
            elif line.startswith("TER"):
                if extracted_lines and extracted_lines[-1][21] == chain_id:
                    extracted_lines.append(line)
            elif line.startswith("END"):
                extracted_lines.append(line)
                break

        with open(output_path, "w") as f:
            f.writelines(extracted_lines)

        print(f"  Extracted chain {chain_id} to: {output_path}")
        return True

    except Exception as e:
        print(f"  ERROR: Failed to extract chain: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Setup target structures")
    parser.add_argument("--output-dir", type=str, default="data/targets", help="Output directory")
    args = parser.parse_args()

    print("=" * 60)
    print("CD3 Binder Pipeline - Setup Target Structures")
    print("=" * 60)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Target structures to download
    targets = [
        {
            "pdb_id": "1XIW",
            "description": "CD3εδ heterodimer",
            "cd3e_chain": "E",
        },
        {
            "pdb_id": "1SY6",
            "description": "CD3εγ + OKT3 Fab complex",
            # NOTE: 1SY6 uses chain A for CD3ε (not E). Chain assignments:
            # - Chain A: CD3ε (epsilon) - the target antigen
            # - Chain H: OKT3 VH (heavy chain variable region)
            # - Chain L: OKT3 VL (light chain variable region)
            "cd3e_chain": "A",
        },
    ]

    success = True

    for target in targets:
        pdb_id = target["pdb_id"]
        print(f"\n{pdb_id}: {target['description']}")

        # Download full PDB
        pdb_path = output_dir / f"{pdb_id.lower()}.pdb"
        if pdb_path.exists():
            print(f"  Already exists: {pdb_path}")
        else:
            if not download_pdb(pdb_id, pdb_path):
                success = False
                continue

        # Extract CD3ε chain
        cd3e_path = output_dir / f"cd3_epsilon_{pdb_id.lower()}.pdb"
        if not extract_cd3_epsilon(pdb_path, cd3e_path, target["cd3e_chain"]):
            success = False

    # Create combined target file for design
    print("\nCreating combined target file...")
    combined_path = output_dir / "cd3_epsilon_combined.pdb"

    # Use 1XIW as primary (simpler structure)
    import shutil
    primary = output_dir / "cd3_epsilon_1xiw.pdb"
    if primary.exists():
        shutil.copy(primary, combined_path)
        print(f"  Primary target: {combined_path}")

    # Extract CD3ε sequence
    print("\nExtracting CD3ε sequence...")

    aa_3to1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

    sequences = {}
    for target in targets:
        pdb_id = target["pdb_id"].lower()
        cd3e_path = output_dir / f"cd3_epsilon_{pdb_id}.pdb"

        if cd3e_path.exists():
            sequence = []
            seen = set()
            with open(cd3e_path, "r") as f:
                for line in f:
                    if line.startswith("ATOM"):
                        res_name = line[17:20].strip()
                        res_num = int(line[22:26])
                        key = (res_num, res_name)
                        if key not in seen and res_name in aa_3to1:
                            seen.add(key)
                            sequence.append((res_num, aa_3to1[res_name]))

            sequence.sort(key=lambda x: x[0])
            seq_str = "".join([aa for _, aa in sequence])
            sequences[pdb_id] = seq_str
            print(f"  {pdb_id}: {len(seq_str)} residues")

    # Save sequences to FASTA
    fasta_path = output_dir / "cd3_epsilon_sequences.fasta"
    with open(fasta_path, "w") as f:
        for pdb_id, seq in sequences.items():
            f.write(f">CD3_epsilon_{pdb_id}\n")
            # Write sequence in 60-character lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    print(f"  Sequences saved to: {fasta_path}")

    print("\n" + "=" * 60)
    if success:
        print("Setup complete!")
        print(f"Target structures ready in: {output_dir}")
    else:
        print("Setup completed with errors. Check output above.")
    print("=" * 60)

    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
