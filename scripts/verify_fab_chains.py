#!/usr/bin/env python3
"""Verify chain IDs in Fab scaffold PDB structures.

Downloads CIF files and identifies actual VH/VL chain labels by parsing
the _struct_asym and _entity categories to map auth chain IDs to CIF chain IDs.

IMPORTANT: RCSB API returns pdbx_strand_id (author IDs like H/L), but BoltzGen
uses struct_asym.id (CIF chain IDs like A/B). These can differ!
"""

import urllib.request
import re
from pathlib import Path


# PDB IDs for each scaffold
SCAFFOLDS = {
    "adalimumab": "6cr1",
    "belimumab": "5y9k",
    "crenezumab": "5vzy",
    "dupilumab": "6wgb",
    "golimumab": "5yoy",
    "guselkumab": "4m6m",
    "mab1": "3h42",
    "necitumumab": "6b3s",
    "nirsevimab": "5udc",
    "sarilumab": "8iow",
    "secukinumab": "6wio",
    "tezepelumab": "5j13",
    "tralokinumab": "5l6y",
    "ustekinumab": "3hmw",
}


def download_cif(pdb_id: str, output_dir: Path) -> Path:
    """Download CIF file from RCSB PDB."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    cif_path = output_dir / f"{pdb_id}.cif"

    if cif_path.exists():
        return cif_path

    print(f"  Downloading {pdb_id.upper()}.cif...")
    with urllib.request.urlopen(url, timeout=30) as response:
        content = response.read()
        cif_path.write_bytes(content)

    return cif_path


def parse_cif_chains(cif_path: Path) -> dict:
    """Parse CIF file to extract chain mappings.

    Returns dict with:
    - struct_asym: list of (cif_chain_id, entity_id) tuples
    - entities: dict of entity_id -> description
    - entity_poly: dict of entity_id -> (pdbx_strand_id, sequence)
    """
    content = cif_path.read_text()

    result = {
        "struct_asym": [],  # (cif_chain_id, entity_id)
        "entities": {},     # entity_id -> description
        "entity_poly": {},  # entity_id -> (auth_chain_id, seq_len)
    }

    # Parse _struct_asym to get CIF chain IDs and their entity mapping
    # Format: _struct_asym.id, _struct_asym.entity_id
    asym_match = re.search(
        r"_struct_asym\.id\s*\n_struct_asym\.pdbx_blank_PDB_chainid_flag\s*\n_struct_asym\.pdbx_modified\s*\n_struct_asym\.entity_id\s*\n_struct_asym\.details\s*\n(.*?)(?=\n#|\nloop_|\n_)",
        content,
        re.DOTALL
    )
    if asym_match:
        lines = asym_match.group(1).strip().split("\n")
        for line in lines:
            parts = line.split()
            if len(parts) >= 4:
                cif_chain_id = parts[0]
                entity_id = parts[3]
                result["struct_asym"].append((cif_chain_id, entity_id))

    # Parse _entity.pdbx_description to get entity descriptions
    entity_match = re.search(
        r"_entity\.id\s*\n_entity\.type\s*\n_entity\.src_method\s*\n_entity\.pdbx_description\s*\n(.*?)(?=\n#|\nloop_|\n_)",
        content,
        re.DOTALL
    )
    if entity_match:
        # Parse multiline entries
        text = entity_match.group(1)
        # Split by entity records
        entries = re.findall(r"(\d+)\s+'?([^']*)'?\s+'?([^']*)'?\s+'([^']*)'", text, re.DOTALL)
        if not entries:
            entries = re.findall(r"(\d+)\s+(\S+)\s+(\S+)\s+;([^;]*);", text, re.DOTALL)
        if not entries:
            # Simpler pattern for single-line descriptions
            for line in text.strip().split("\n"):
                if line.strip():
                    parts = line.split()
                    if len(parts) >= 4:
                        eid = parts[0]
                        desc = " ".join(parts[3:]).strip("'\"")
                        result["entities"][eid] = desc
        else:
            for entry in entries:
                eid = entry[0]
                desc = entry[3] if len(entry) > 3 else entry[1]
                result["entities"][eid] = desc.strip()

    # Parse _entity_poly to get auth chain IDs (pdbx_strand_id)
    poly_match = re.search(
        r"_entity_poly\.pdbx_strand_id\s*\n_entity_poly\.pdbx_target_identifier\s*\n(.*?)(?=\n#|\nloop_|\n_[a-z])",
        content,
        re.DOTALL
    )
    if poly_match:
        # This is trickier - need to match entity records
        pass

    # Alternative: find pdbx_strand_id values after sequence blocks
    strand_matches = re.findall(r";\s*\n([A-Z0-9,]+)\s+\?", content)
    entity_id = 1
    for strand_id in strand_matches:
        result["entity_poly"][str(entity_id)] = strand_id.split(",")[0]
        entity_id += 1

    return result


def identify_vh_vl_from_cif(cif_path: Path) -> dict:
    """Identify VH and VL chains from CIF file.

    Returns dict with vh_chain and vl_chain (CIF chain IDs).
    """
    result = {"vh_chain": None, "vl_chain": None, "chains": []}

    chain_info = parse_cif_chains(cif_path)

    # Try to match entity descriptions to heavy/light chains
    heavy_entity = None
    light_entity = None

    for entity_id, desc in chain_info["entities"].items():
        desc_lower = desc.lower()
        if "heavy" in desc_lower or "vh " in desc_lower:
            heavy_entity = entity_id
        elif "light" in desc_lower or "vl " in desc_lower or "kappa" in desc_lower or "lambda" in desc_lower:
            light_entity = entity_id

    # Map entity to CIF chain ID
    for cif_chain, entity_id in chain_info["struct_asym"]:
        entity_desc = chain_info["entities"].get(entity_id, "")
        auth_chain = chain_info["entity_poly"].get(entity_id, "?")

        result["chains"].append({
            "cif_chain": cif_chain,
            "entity_id": entity_id,
            "auth_chain": auth_chain,
            "description": entity_desc[:60],
        })

        if entity_id == heavy_entity and result["vh_chain"] is None:
            result["vh_chain"] = cif_chain
        elif entity_id == light_entity and result["vl_chain"] is None:
            result["vl_chain"] = cif_chain

    return result


def main():
    print("=" * 70)
    print("Verifying Fab Scaffold Chain IDs (CIF struct_asym.id)")
    print("=" * 70)
    print()
    print("NOTE: BoltzGen uses struct_asym.id (e.g., A, B) NOT pdbx_strand_id (e.g., H, L)")
    print()

    # Download CIF files to temp dir
    tmp_dir = Path("/tmp/fab_scaffold_verify")
    tmp_dir.mkdir(exist_ok=True)

    results = {}

    for name, pdb_id in SCAFFOLDS.items():
        print(f"Checking {name} ({pdb_id.upper()})...")
        try:
            cif_path = download_cif(pdb_id, tmp_dir)
            info = identify_vh_vl_from_cif(cif_path)
            results[name] = info

            print(f"  VH chain (CIF): {info['vh_chain']}")
            print(f"  VL chain (CIF): {info['vl_chain']}")
            for c in info["chains"][:4]:
                print(f"    {c['cif_chain']}: entity={c['entity_id']}, auth={c['auth_chain']}, {c['description']}")

        except Exception as e:
            print(f"  Error: {e}")
            results[name] = {"vh_chain": None, "vl_chain": None, "error": str(e)}
        print()

    # Print summary
    print("=" * 70)
    print("SUMMARY - CIF chain IDs for setup_fab_scaffolds.py")
    print("=" * 70)
    print(f"{'Scaffold':<15} {'PDB':<6} {'VH':<4} {'VL':<4}")
    print("-" * 35)

    for name, pdb_id in SCAFFOLDS.items():
        info = results[name]
        vh = info.get("vh_chain") or "?"
        vl = info.get("vl_chain") or "?"
        print(f"{name:<15} {pdb_id.upper():<6} {vh:<4} {vl:<4}")

    print()
    print("Python dict for FAB_SCAFFOLDS:")
    print("-" * 70)
    for name, pdb_id in SCAFFOLDS.items():
        info = results[name]
        vh = info.get("vh_chain")
        vl = info.get("vl_chain")
        if vh and vl:
            print(f'    "{name}": {{"vh_chain": "{vh}", "vl_chain": "{vl}"}},')
        else:
            print(f'    # "{name}": NEEDS MANUAL CHECK')

    return 0


if __name__ == "__main__":
    exit(main())
