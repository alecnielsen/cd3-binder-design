#!/usr/bin/env python3
"""Setup Fab scaffold files for BoltzGen CDR redesign.

Downloads CIF files from RCSB PDB and creates YAML scaffold specifications
with CDR region definitions for BoltzGen Fab CDR redesign.

Usage:
    python scripts/setup_fab_scaffolds.py [--output-dir data/fab_scaffolds]
"""

import argparse
import urllib.request
from pathlib import Path


# Fab scaffold definitions with CDR regions (Kabat numbering)
# Each scaffold has:
#   - pdb_id: PDB ID to download
#   - vh_chain: Chain ID for VH
#   - vl_chain: Chain ID for VL
#   - vh_range: Residue range for VH (1-indexed)
#   - vl_range: Residue range for VL (1-indexed)
#   - cdr_h1/h2/h3: H-CDR residue ranges
#   - cdr_l1/l2/l3: L-CDR residue ranges
FAB_SCAFFOLDS = {
    "adalimumab": {
        "pdb_id": "6cr1",
        "vh_chain": "B",
        "vl_chain": "A",
        "vh_range": (1, 121),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 110),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "belimumab": {
        "pdb_id": "5y9k",
        "vh_chain": "B",
        "vl_chain": "A",
        "vh_range": (1, 120),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 108),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "crenezumab": {
        "pdb_id": "5vzy",
        "vh_chain": "A",
        "vl_chain": "B",
        "vh_range": (1, 121),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 110),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "dupilumab": {
        "pdb_id": "6wgb",
        "vh_chain": "A",
        "vl_chain": "B",
        "vh_range": (1, 120),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 109),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "golimumab": {
        "pdb_id": "5yoy",
        "vh_chain": "G",
        "vl_chain": "D",
        "vh_range": (1, 121),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 110),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "guselkumab": {
        "pdb_id": "4m6m",
        "vh_chain": "B",
        "vl_chain": "A",
        "vh_range": (1, 121),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 110),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "mab1": {
        "pdb_id": "3h42",
        "vh_chain": "D",
        "vl_chain": "C",
        "vh_range": (1, 119),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 107),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "necitumumab": {
        "pdb_id": "6b3s",
        "vh_chain": "B",
        "vl_chain": "C",
        "vh_range": (1, 118),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 106),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "nirsevimab": {
        "pdb_id": "5udc",
        "vh_chain": "A",
        "vl_chain": "B",
        "vh_range": (1, 121),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 110),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "sarilumab": {
        "pdb_id": "8iow",
        "vh_chain": "D",
        "vl_chain": "C",
        "vh_range": (1, 120),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 109),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "secukinumab": {
        "pdb_id": "6wio",
        "vh_chain": "A",
        "vl_chain": "B",
        "vh_range": (1, 120),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 108),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "tezepelumab": {
        "pdb_id": "5j13",
        "vh_chain": "C",
        "vl_chain": "B",
        "vh_range": (1, 121),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 109),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "tralokinumab": {
        "pdb_id": "5l6y",
        "vh_chain": "B",
        "vl_chain": "C",
        "vh_range": (1, 120),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 109),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
    "ustekinumab": {
        "pdb_id": "3hmw",
        "vh_chain": "B",
        "vl_chain": "A",
        "vh_range": (1, 119),
        "vl_range": (1, 107),
        "cdr_h1": (26, 32),
        "cdr_h2": (52, 57),
        "cdr_h3": (99, 108),
        "cdr_l1": (24, 34),
        "cdr_l2": (50, 56),
        "cdr_l3": (89, 97),
    },
}


def download_cif(pdb_id: str, output_path: Path) -> bool:
    """Download CIF file from RCSB PDB.

    Args:
        pdb_id: 4-character PDB ID (lowercase).
        output_path: Path to write CIF file.

    Returns:
        True if successful, False otherwise.
    """
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"

    try:
        print(f"  Downloading {pdb_id.upper()}.cif from RCSB...")
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read()
            output_path.write_bytes(content)
        return True
    except Exception as e:
        print(f"  ERROR: Failed to download {pdb_id}: {e}")
        return False


def create_scaffold_yaml(name: str, scaffold: dict, output_dir: Path) -> bool:
    """Create BoltzGen scaffold YAML file.

    Args:
        name: Scaffold name (e.g., "adalimumab").
        scaffold: Scaffold definition dict.
        output_dir: Directory to write YAML file.

    Returns:
        True if successful, False otherwise.
    """
    pdb_id = scaffold["pdb_id"]
    vh_chain = scaffold["vh_chain"]
    vl_chain = scaffold["vl_chain"]
    vh_start, vh_end = scaffold["vh_range"]
    vl_start, vl_end = scaffold["vl_range"]

    # CDR regions
    h1_start, h1_end = scaffold["cdr_h1"]
    h2_start, h2_end = scaffold["cdr_h2"]
    h3_start, h3_end = scaffold["cdr_h3"]
    l1_start, l1_end = scaffold["cdr_l1"]
    l2_start, l2_end = scaffold["cdr_l2"]
    l3_start, l3_end = scaffold["cdr_l3"]

    # Build YAML content
    # Using exact BoltzGen format from examples
    yaml_content = f"""# BoltzGen Fab scaffold: {name}
# PDB: {pdb_id.upper()}
# Human antibody framework with CDR redesign regions

path: {name}.{pdb_id}.cif

include:
  - chain:
      id: {vh_chain}
      res_index: {vh_start}..{vh_end}
  - chain:
      id: {vl_chain}
      res_index: {vl_start}..{vl_end}

design:
  - chain:
      id: {vh_chain}
      res_index: {h1_start}..{h1_end},{h2_start}..{h2_end},{h3_start}..{h3_end}
  - chain:
      id: {vl_chain}
      res_index: {l1_start}..{l1_end},{l2_start}..{l2_end},{l3_start}..{l3_end}

design_insertions:
  - insertion:
      id: {vh_chain}
      res_index: {h3_start}
      num_residues: 3..21
  - insertion:
      id: {vl_chain}
      res_index: {l3_start}
      num_residues: 3..12
"""

    yaml_path = output_dir / f"{name}.{pdb_id}.yaml"
    yaml_path.write_text(yaml_content)
    return True


def main():
    parser = argparse.ArgumentParser(description="Setup Fab scaffolds for BoltzGen")
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data/fab_scaffolds",
        help="Output directory for scaffold files",
    )
    parser.add_argument(
        "--scaffolds",
        type=str,
        nargs="+",
        default=None,
        help="Specific scaffolds to download (default: all)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    scaffolds = args.scaffolds or list(FAB_SCAFFOLDS.keys())

    print("=" * 60)
    print("BoltzGen Fab Scaffold Setup")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    print(f"Scaffolds: {', '.join(scaffolds)}")
    print()

    success_count = 0
    failed = []

    for name in scaffolds:
        if name not in FAB_SCAFFOLDS:
            print(f"WARNING: Unknown scaffold '{name}', skipping")
            failed.append(name)
            continue

        scaffold = FAB_SCAFFOLDS[name]
        pdb_id = scaffold["pdb_id"]
        cif_path = output_dir / f"{name}.{pdb_id}.cif"

        print(f"Processing {name} ({pdb_id.upper()})...")

        # Download CIF if needed
        if cif_path.exists():
            print(f"  CIF already exists: {cif_path}")
        else:
            if not download_cif(pdb_id, cif_path):
                failed.append(name)
                continue

        # Create YAML scaffold spec
        if create_scaffold_yaml(name, scaffold, output_dir):
            print(f"  Created YAML: {name}.{pdb_id}.yaml")
            success_count += 1
        else:
            failed.append(name)

    print()
    print("=" * 60)
    print(f"Setup complete: {success_count}/{len(scaffolds)} scaffolds")
    if failed:
        print(f"Failed: {', '.join(failed)}")
    print("=" * 60)

    # List files
    print(f"\nFiles in {output_dir}:")
    for f in sorted(output_dir.glob("*")):
        print(f"  {f.name}")

    return 0 if not failed else 1


if __name__ == "__main__":
    exit(main())
