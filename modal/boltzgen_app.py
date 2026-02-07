"""BoltzGen Modal deployment for de novo binder design.

This module provides a Modal app for running BoltzGen
on GPU infrastructure for de novo antibody design.

IMPORTANT: This module extracts ONLY the target chain from input PDBs
before passing to BoltzGen. Multi-chain PDBs (e.g., with bound antibodies)
would bias design if passed in full.

Deploy with:
    modal deploy modal/boltzgen_app.py

Run locally:
    # First download model weights
    modal run modal/boltzgen_app.py --download

    # Then run design
    modal run modal/boltzgen_app.py --target-pdb data/targets/cd3.pdb --target-chain A
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import modal

MINUTES = 60

app = modal.App("boltzgen-cd3")

# Container with BoltzGen - install from GitHub main for latest fixes
# (PyPI 0.2.0 is from Dec 10, but refolding fix landed Dec 17)
boltzgen_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git")
    .pip_install("boltzgen @ git+https://github.com/HannesStark/boltzgen.git@main")
)

# Persistent volume for model weights
boltzgen_model_volume = modal.Volume.from_name("boltzgen-models", create_if_missing=True)
models_dir = Path("/models/boltzgen")

# Image for downloading model
download_image = (
    modal.Image.debian_slim()
    .pip_install("huggingface-hub==0.36.0", "hf_transfer")
    .env({"HF_HUB_ENABLE_HF_TRANSFER": "1"})
)


def extract_chain_from_pdb_content(pdb_content: str, chain_id: str) -> str:
    """Extract a single chain from PDB content.

    CRITICAL: This function ensures only the target chain is passed to BoltzGen.
    Passing multi-chain PDBs (e.g., CD3 + bound antibody) can bias designs toward
    the pre-occupied epitope or generate antibody-like sequences.

    Args:
        pdb_content: Full PDB file content as string.
        chain_id: Chain ID to extract (e.g., "A" for CD3ε).

    Returns:
        PDB string containing only the specified chain.

    Raises:
        ValueError: If the specified chain is not found or has no atoms.
    """
    lines = []
    atom_count = 0

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if len(line) > 21 and line[21] == chain_id:
                lines.append(line)
                atom_count += 1
        elif line.startswith("TER"):
            if lines and len(lines[-1]) > 21 and lines[-1][21] == chain_id:
                lines.append(line)
        elif line.startswith("END"):
            lines.append(line)
            break

    if atom_count == 0:
        available_chains = set()
        for line in pdb_content.split("\n"):
            if line.startswith("ATOM") and len(line) > 21:
                available_chains.add(line[21])
        raise ValueError(
            f"Chain '{chain_id}' not found in PDB. "
            f"Available chains: {sorted(available_chains)}"
        )

    return "\n".join(lines)


def parse_pdb_sequence(pdb_content: str, chain_id: str = "A") -> str:
    """Extract amino acid sequence from PDB content."""
    aa_3to1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

    sequence = []
    seen_residues = set()

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM"):
            chain = line[21]
            if chain != chain_id:
                continue

            res_name = line[17:20].strip()
            res_num = int(line[22:26])
            insertion_code = line[26] if len(line) > 26 else " "
            res_key = (res_num, insertion_code, res_name)

            if res_key not in seen_residues and res_name in aa_3to1:
                seen_residues.add(res_key)
                sequence.append((res_num, insertion_code, aa_3to1[res_name]))

    sequence.sort(key=lambda x: (x[0], x[1]))
    return "".join([aa for _, _, aa in sequence])


def build_design_spec_yaml(
    target_cif_path: str,
    target_chain: str = "A",
    binder_length: int | tuple[int, int] = 120,
    hotspot_residues: list[int] | None = None,
) -> str:
    """Build BoltzGen design specification YAML.

    Uses the correct BoltzGen format:
    - Target referenced via file path (not inline sequence)
    - Binder specified with range notation (e.g., "80..140")
    - Binding sites via binding_types section

    Args:
        target_cif_path: Path to target CIF file (relative to YAML location).
        target_chain: Chain ID to use from target file.
        binder_length: Length of binder - either int or (min, max) tuple.
        hotspot_residues: Optional list of target residues to focus binding on.

    Returns:
        YAML string for BoltzGen design spec.
    """
    # Format binder length as range notation
    if isinstance(binder_length, tuple):
        min_len, max_len = binder_length
        length_spec = f"{min_len}..{max_len}"
    else:
        # Use a small range around the target length for flexibility
        length_spec = f"{binder_length - 10}..{binder_length + 10}"

    # Build YAML manually (avoid PyYAML dependency)
    lines = ["entities:"]

    # Designed binder chain - use range notation
    lines.append("  - protein:")
    lines.append("      id: B")
    lines.append(f"      sequence: {length_spec}")

    # Target from file
    lines.append("  - file:")
    lines.append(f"      path: {target_cif_path}")
    lines.append("      include:")
    lines.append("        - chain:")
    lines.append(f"            id: {target_chain}")

    # Add binding site specification if hotspot residues provided
    if hotspot_residues:
        lines.append("")
        lines.append("binding_types:")
        lines.append("  - chain:")
        lines.append(f"      id: {target_chain}")
        residues_str = ",".join(str(r) for r in hotspot_residues)
        lines.append(f"      binding: {residues_str}")

    return "\n".join(lines)


@app.function(
    volumes={models_dir: boltzgen_model_volume},
    timeout=60 * MINUTES,
    image=download_image,
)
def download_model(force_download: bool = False):
    """Download BoltzGen model weights to Modal volume."""
    from huggingface_hub import snapshot_download

    # All BoltzGen checkpoints are in boltzgen/boltzgen-1
    # See: https://huggingface.co/boltzgen/boltzgen-1
    print("Downloading BoltzGen model weights from boltzgen/boltzgen-1...")
    snapshot_download(
        repo_id="boltzgen/boltzgen-1",
        local_dir=models_dir / "boltzgen-1",
        force_download=force_download,
    )

    # Also download inference data (molecular info) - this is a dataset repo
    print("Downloading inference data...")
    snapshot_download(
        repo_id="boltzgen/inference-data",
        repo_type="dataset",
        local_dir=models_dir / "inference-data",
        force_download=force_download,
    )

    boltzgen_model_volume.commit()
    print(f"Models downloaded to {models_dir}")


@app.function(
    image=boltzgen_image,
    volumes={models_dir: boltzgen_model_volume},
    timeout=60 * MINUTES,
    gpu="A100",
)
def run_boltzgen(
    target_pdb_content: str,
    target_chain: str = "A",
    num_designs: int = 100,
    binder_length: int = 120,
    hotspot_residues: list[int] | None = None,
    protocol: str = "nanobody-anything",
    seed: int = 42,
) -> list[dict]:
    """Run BoltzGen to design binders.

    IMPORTANT: This function extracts ONLY the target chain before passing
    to BoltzGen. Multi-chain PDBs (e.g., CD3 + bound Fab) would bias designs
    if passed in full.

    Args:
        target_pdb_content: PDB file content as string (can be multi-chain).
        target_chain: Chain ID of target in PDB (will be extracted).
        num_designs: Number of designs to generate.
        binder_length: Length of binder (120 for VHH, ~250 for scFv).
        hotspot_residues: Optional residues to target for binding.
        protocol: BoltzGen protocol (nanobody-anything, protein-anything).
        seed: Random seed for reproducibility.

    Returns:
        List of design dictionaries with sequences and scores.

    Raises:
        ValueError: If target_chain is not found in the PDB.
    """
    # CRITICAL: Extract only the target chain to avoid design bias
    extracted_pdb = extract_chain_from_pdb_content(target_pdb_content, target_chain)

    # Get target sequence for logging
    target_sequence = parse_pdb_sequence(extracted_pdb, target_chain)
    print(f"Target sequence ({len(target_sequence)} aa): {target_sequence[:50]}...")

    # Write target PDB to file (BoltzGen needs file reference, not inline sequence)
    # Note: BoltzGen accepts both PDB and CIF files
    target_pdb_path = Path("target.pdb")
    target_pdb_path.write_text(extracted_pdb)
    print(f"Wrote extracted target chain {target_chain} to {target_pdb_path}")

    # Build design spec YAML with file reference
    spec_yaml = build_design_spec_yaml(
        target_cif_path="target.pdb",  # BoltzGen accepts PDB too
        target_chain=target_chain,
        binder_length=binder_length,
        hotspot_residues=hotspot_residues,
    )

    spec_path = Path("design_spec.yaml")
    spec_path.write_text(spec_yaml)
    print(f"Design spec:\n{spec_yaml}")

    output_dir = Path("boltzgen_output")
    output_dir.mkdir(exist_ok=True)

    # Run BoltzGen CLI
    cmd = [
        "boltzgen", "run", str(spec_path),
        "--protocol", protocol,
        "--output", str(output_dir),
        "--num_designs", str(num_designs),
        "--cache", str(models_dir),
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    print(f"STDOUT: {result.stdout}")
    if result.returncode != 0:
        print(f"STDERR: {result.stderr}")
        # Fallback: check if intermediate designs exist despite failure
        # (kept for robustness even though Dec 17 refolding fix should resolve this)
        ifold_dir = output_dir / "intermediate_designs_inverse_folded"
        if ifold_dir.exists():
            print(f"Pipeline failed but intermediate designs found at {ifold_dir}")
        else:
            raise RuntimeError(f"BoltzGen failed: {result.stderr}")

    # Parse outputs - BoltzGen outputs FASTA and/or CIF files
    designs = []

    # Look for output files - prefer final_ranked_designs, fall back to intermediate
    search_dirs = [
        output_dir / "final_ranked_designs",
        output_dir / "intermediate_designs_inverse_folded",
        output_dir / "intermediate_designs",
    ]

    fasta_files = []
    cif_files = []
    json_files = []
    csv_files = []

    for search_dir in search_dirs:
        if search_dir.exists():
            fasta_files.extend(list(search_dir.glob("**/*.fasta")))
            fasta_files.extend(list(search_dir.glob("**/*.fa")))
            cif_files.extend(list(search_dir.glob("**/*.cif")))
            json_files.extend(list(search_dir.glob("**/*.json")))
            csv_files.extend(list(search_dir.glob("**/*.csv")))

    print(f"Found {len(fasta_files)} FASTA, {len(cif_files)} CIF, {len(json_files)} JSON, {len(csv_files)} CSV files")

    # Parse metrics from CSV if available (prefer metrics.csv which has design scores)
    # IMPORTANT: CSV row order = BoltzGen's internal ranking (decision tree output)
    metrics_by_design = {}
    # Sort to prefer metrics.csv over other CSVs
    csv_files_sorted = sorted(csv_files, key=lambda f: "metrics" in f.name, reverse=True)
    for csv_file in csv_files_sorted:
        try:
            content = csv_file.read_text()
            csv_lines = content.strip().split("\n")
            if len(csv_lines) < 2:
                continue
            headers = csv_lines[0].split(",")
            print(f"  CSV {csv_file.name}: {headers[:8]}...")
            for rank_idx, row in enumerate(csv_lines[1:]):
                values = row.split(",")
                if len(values) >= len(headers):
                    row_dict = dict(zip(headers, values))
                    # Use 'id' column as key (BoltzGen uses this)
                    design_id = row_dict.get("id", row_dict.get("file_name", ""))
                    if design_id and design_id not in metrics_by_design:
                        # Row order in CSV = BoltzGen's rank (0-indexed)
                        row_dict["_boltzgen_rank"] = rank_idx + 1
                        metrics_by_design[design_id] = row_dict
                        iptm = row_dict.get("design_to_target_iptm", "N/A")
                        ptm = row_dict.get("design_ptm", "N/A")
                        print(f"    rank {rank_idx + 1}: {design_id}: ipTM={iptm}, pTM={ptm}")
        except Exception as e:
            print(f"  Error parsing CSV {csv_file}: {e}")

    # Parse FASTA outputs for sequences
    for fasta_file in fasta_files:
        content = fasta_file.read_text()
        # Simple FASTA parser
        current_header = None
        current_seq = []

        for line in content.split("\n"):
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_seq:
                    designs.append({
                        "sequence": "".join(current_seq),
                        "header": current_header,
                        "source_file": fasta_file.name,
                    })
                current_header = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line)

        # Don't forget last sequence
        if current_header and current_seq:
            designs.append({
                "sequence": "".join(current_seq),
                "header": current_header,
                "source_file": fasta_file.name,
            })

    # If no FASTA files, try to extract sequences from CIF files
    if not designs and cif_files:
        print("No FASTA files found, extracting sequences from CIF files...")
        # Only process inverse-folded files (skip refold_cif to avoid duplicates)
        ifold_cifs = [
            f for f in cif_files
            if "inverse_folded" in str(f.parent) and "refold_cif" not in str(f)
        ]
        if not ifold_cifs:
            ifold_cifs = cif_files  # Fallback to all CIF files
        print(f"  Processing {len(ifold_cifs)} CIF files...")
        for cif_file in ifold_cifs:
            # Skip the design_spec.cif visualization file (only target, no design)
            if cif_file.name == "design_spec.cif":
                continue
            print(f"  Processing {cif_file}...")
            try:
                content = cif_file.read_text()
                lines = content.split("\n")

                # Parse _entity_poly loop for sequences
                # Format: entity_id type strand_id sequence
                # Entity 1 is the binder (chain B), Entity 2 is the target
                in_entity_poly = False
                entity_poly_cols = []
                seq_col_idx = -1
                binder_sequence = None

                for line in lines:
                    line = line.strip()
                    if line.startswith("_entity_poly."):
                        in_entity_poly = True
                        col_name = line.split(".")[1]
                        entity_poly_cols.append(col_name)
                        if col_name == "pdbx_seq_one_letter_code":
                            seq_col_idx = len(entity_poly_cols) - 1
                    elif in_entity_poly:
                        if line.startswith("_") or line.startswith("#") or line == "":
                            in_entity_poly = False
                            continue
                        if line.startswith("loop_"):
                            continue
                        # Data line - parse entity_id and sequence
                        parts = line.split()
                        if len(parts) >= 4 and seq_col_idx >= 0:
                            entity_id = parts[0]
                            # Sequence is typically the last column
                            sequence = parts[-1] if len(parts) > seq_col_idx else parts[seq_col_idx]
                            print(f"    Entity {entity_id}: {sequence[:50]}... ({len(sequence)} aa)")
                            # Entity 1 is the binder (first in YAML = chain B)
                            if entity_id == "1" and "X" not in sequence:
                                binder_sequence = sequence

                if binder_sequence:
                    designs.append({
                        "sequence": binder_sequence,
                        "header": cif_file.stem,
                        "source_file": cif_file.name,
                    })
                    print(f"    -> Extracted binder: {binder_sequence[:40]}...")
            except Exception as e:
                print(f"    Error parsing {cif_file.name}: {e}")
                import traceback
                traceback.print_exc()

    # Try to parse confidence scores from JSON if available
    for json_file in json_files:
        try:
            scores = json.loads(json_file.read_text())
            # Match scores to designs by index if possible
            if isinstance(scores, list):
                for i, score in enumerate(scores):
                    if i < len(designs):
                        designs[i].update({
                            "confidence": score.get("confidence", 0.0),
                            "plddt": score.get("plddt", 0.0),
                            "ptm": score.get("ptm", 0.0),
                        })
        except (json.JSONDecodeError, KeyError):
            pass

    # Add design indices and merge metrics
    for i, design in enumerate(designs):
        design["design_idx"] = i
        # Try to match with CSV metrics by design name
        design_name = design.get("header", "")
        if design_name in metrics_by_design:
            metrics = metrics_by_design[design_name]
            # BoltzGen uses these column names
            try:
                design["ipTM"] = float(metrics.get("design_to_target_iptm", 0))
                design["pTM"] = float(metrics.get("design_ptm", 0))
                design["pae_min"] = float(metrics.get("min_design_to_target_pae", 0))
                design["rmsd"] = float(metrics.get("filter_rmsd", 0))
                design["rmsd_design"] = float(metrics.get("filter_rmsd_design", 0))
                design["boltzgen_rank"] = metrics.get("_boltzgen_rank", i + 1)
            except (ValueError, TypeError):
                pass  # Skip if conversion fails
        else:
            # Fallback: use parse order as rank
            design["boltzgen_rank"] = i + 1

    print(f"Parsed {len(designs)} designs")
    return designs


# Available Fab scaffolds for CDR redesign
# Each scaffold must have corresponding CIF and YAML files in data/fab_scaffolds/
AVAILABLE_FAB_SCAFFOLDS = [
    "adalimumab",
    "belimumab",
    "crenezumab",
    "dupilumab",
    "golimumab",
    "guselkumab",
    "mab1",
    "necitumumab",
    "nirsevimab",
    "sarilumab",
    "secukinumab",
    "tezepelumab",
    "tralokinumab",
    "ustekinumab",
]


def build_fab_target_yaml(
    target_cif_path: str,
    target_chain: str,
    scaffold_yaml_paths: list[str],
) -> str:
    """Build BoltzGen target YAML for Fab CDR redesign.

    This creates a YAML spec that references:
    - The target structure (CD3)
    - Multiple scaffold YAML files (each defines a Fab with CDR design regions)

    Args:
        target_cif_path: Path to target CIF/PDB file.
        target_chain: Chain ID of target to use.
        scaffold_yaml_paths: List of paths to scaffold YAML files.

    Returns:
        YAML string for BoltzGen target spec.
    """
    lines = ["entities:"]

    # Target from file
    lines.append("  - file:")
    lines.append(f"      path: {target_cif_path}")
    lines.append("      include:")
    lines.append("        - chain:")
    lines.append(f"            id: {target_chain}")

    # Scaffold files (each defines a Fab)
    lines.append("  - file:")
    lines.append("      path:")
    for scaffold_path in scaffold_yaml_paths:
        lines.append(f"        - {scaffold_path}")

    return "\n".join(lines)


@app.function(
    image=boltzgen_image,
    volumes={models_dir: boltzgen_model_volume},
    timeout=90 * MINUTES,
    gpu="A100",
)
def run_boltzgen_fab(
    target_pdb_content: str,
    target_chain: str = "A",
    scaffold_names: list[str] | None = None,
    scaffold_yaml_contents: dict[str, str] | None = None,
    scaffold_cif_contents: dict[str, str] | None = None,
    num_designs: int = 100,
    seed: int = 42,
) -> list[dict]:
    """Run BoltzGen Fab CDR redesign.

    Designs Fab binders by redesigning CDR loops while maintaining
    human framework regions from proven therapeutic antibodies.

    Args:
        target_pdb_content: PDB file content as string (can be multi-chain).
        target_chain: Chain ID of target in PDB (will be extracted).
        scaffold_names: Names of scaffolds to use (e.g., ["adalimumab"]).
        scaffold_yaml_contents: Dict mapping scaffold names to YAML content.
        scaffold_cif_contents: Dict mapping scaffold names to CIF content.
        num_designs: Number of designs to generate per scaffold.
        seed: Random seed for reproducibility.

    Returns:
        List of design dictionaries with VH/VL sequences and scores.

    Raises:
        ValueError: If target_chain is not found or scaffolds are missing.
    """
    # CRITICAL: Extract only the target chain to avoid design bias
    extracted_pdb = extract_chain_from_pdb_content(target_pdb_content, target_chain)

    # Get target sequence for logging
    target_sequence = parse_pdb_sequence(extracted_pdb, target_chain)
    print(f"Target sequence ({len(target_sequence)} aa): {target_sequence[:50]}...")

    # Write target PDB to file
    target_pdb_path = Path("target.pdb")
    target_pdb_path.write_text(extracted_pdb)
    print(f"Wrote extracted target chain {target_chain} to {target_pdb_path}")

    # Create scaffold directory and write files
    scaffold_dir = Path("fab_scaffolds")
    scaffold_dir.mkdir(exist_ok=True)

    if not scaffold_names:
        scaffold_names = ["adalimumab"]  # Default to adalimumab

    scaffold_yaml_paths = []
    for name in scaffold_names:
        if scaffold_yaml_contents and name in scaffold_yaml_contents:
            yaml_content = scaffold_yaml_contents[name]
        else:
            raise ValueError(f"Missing YAML content for scaffold: {name}")

        if scaffold_cif_contents and name in scaffold_cif_contents:
            cif_content = scaffold_cif_contents[name]
        else:
            raise ValueError(f"Missing CIF content for scaffold: {name}")

        # Find CIF filename from YAML (first line with 'path:')
        cif_filename = None
        for line in yaml_content.split("\n"):
            if line.strip().startswith("path:") and ".cif" in line:
                cif_filename = line.split(":")[-1].strip()
                break

        if not cif_filename:
            cif_filename = f"{name}.cif"

        # Write CIF file
        cif_path = scaffold_dir / cif_filename
        cif_path.write_text(cif_content)
        print(f"Wrote scaffold CIF: {cif_path}")

        # Write YAML file (path stays as cif_filename since YAML and CIF are in same dir)
        yaml_path = scaffold_dir / f"{name}.yaml"
        # Keep relative path - YAML and CIF are in the same directory
        updated_yaml = yaml_content  # No path update needed - already relative
        yaml_path.write_text(updated_yaml)
        print(f"Wrote scaffold YAML: {yaml_path}")

        scaffold_yaml_paths.append(str(yaml_path))

    # Build target spec YAML
    target_yaml = build_fab_target_yaml(
        target_cif_path="target.pdb",
        target_chain=target_chain,
        scaffold_yaml_paths=scaffold_yaml_paths,
    )

    spec_path = Path("fab_design_spec.yaml")
    spec_path.write_text(target_yaml)
    print(f"Design spec:\n{target_yaml}")

    output_dir = Path("boltzgen_fab_output")
    output_dir.mkdir(exist_ok=True)

    # Run BoltzGen CLI with antibody-anything protocol
    cmd = [
        "boltzgen", "run", str(spec_path),
        "--protocol", "antibody-anything",
        "--output", str(output_dir),
        "--num_designs", str(num_designs),
        "--cache", str(models_dir),
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    print(f"STDOUT: {result.stdout}")
    if result.returncode != 0:
        print(f"STDERR: {result.stderr}")
        # Fallback: check if intermediate designs exist despite failure
        # (kept for robustness even though Dec 17 refolding fix should resolve this)
        ifold_dir = output_dir / "intermediate_designs_inverse_folded"
        if ifold_dir.exists():
            print(f"Pipeline failed but intermediate designs found at {ifold_dir}")
        else:
            raise RuntimeError(f"BoltzGen Fab failed: {result.stderr}")

    # Parse outputs - look for Fab designs with VH/VL chains
    designs = []

    # Look for output files
    search_dirs = [
        output_dir / "final_ranked_designs",
        output_dir / "intermediate_designs_inverse_folded",
        output_dir / "intermediate_designs",
    ]

    fasta_files = []
    cif_files = []
    csv_files = []

    for search_dir in search_dirs:
        if search_dir.exists():
            fasta_files.extend(list(search_dir.glob("**/*.fasta")))
            fasta_files.extend(list(search_dir.glob("**/*.fa")))
            cif_files.extend(list(search_dir.glob("**/*.cif")))
            csv_files.extend(list(search_dir.glob("**/*.csv")))

    print(f"Found {len(fasta_files)} FASTA, {len(cif_files)} CIF, {len(csv_files)} CSV files")

    # Parse metrics from CSV
    # IMPORTANT: CSV row order = BoltzGen's internal ranking (decision tree output)
    metrics_by_design = {}
    csv_files_sorted = sorted(csv_files, key=lambda f: "metrics" in f.name, reverse=True)
    for csv_file in csv_files_sorted:
        try:
            content = csv_file.read_text()
            csv_lines = content.strip().split("\n")
            if len(csv_lines) < 2:
                continue
            headers = csv_lines[0].split(",")
            print(f"  CSV {csv_file.name}: {headers[:8]}...")
            for rank_idx, row in enumerate(csv_lines[1:]):
                values = row.split(",")
                if len(values) >= len(headers):
                    row_dict = dict(zip(headers, values))
                    design_id = row_dict.get("id", row_dict.get("file_name", ""))
                    if design_id and design_id not in metrics_by_design:
                        row_dict["_boltzgen_rank"] = rank_idx + 1
                        metrics_by_design[design_id] = row_dict
        except Exception as e:
            print(f"  Error parsing CSV {csv_file}: {e}")

    # Parse CIF files to extract VH/VL sequences
    # For Fab designs, we expect two chains: VH and VL
    ifold_cifs = [
        f for f in cif_files
        if "inverse_folded" in str(f.parent) and "refold_cif" not in str(f)
    ]
    if not ifold_cifs:
        ifold_cifs = cif_files

    print(f"Processing {len(ifold_cifs)} CIF files for VH/VL extraction...")

    for cif_file in ifold_cifs:
        if cif_file.name == "fab_design_spec.cif":
            continue

        print(f"  Processing {cif_file}...")
        try:
            content = cif_file.read_text()
            lines = content.split("\n")

            # Parse _entity_poly for sequences
            # For Fab: Entity 1 = VH, Entity 2 = VL, Entity 3 = target
            in_entity_poly = False
            entity_poly_cols = []
            sequences_by_entity = {}

            for line in lines:
                line_stripped = line.strip()
                if line_stripped.startswith("_entity_poly."):
                    in_entity_poly = True
                    col_name = line_stripped.split(".")[1]
                    entity_poly_cols.append(col_name)
                elif in_entity_poly:
                    if line_stripped.startswith("_") or line_stripped.startswith("#") or line_stripped == "":
                        in_entity_poly = False
                        continue
                    if line_stripped.startswith("loop_"):
                        continue
                    parts = line_stripped.split()
                    if len(parts) >= 2:
                        entity_id = parts[0]
                        # Sequence is typically the last column
                        sequence = parts[-1] if len(parts) > 3 else ""
                        if sequence and "X" not in sequence:
                            sequences_by_entity[entity_id] = sequence
                            print(f"    Entity {entity_id}: {sequence[:40]}... ({len(sequence)} aa)")

            # Identify VH and VL by sequence characteristics
            # VH: starts with EVQL/QVQL/DVQL (heavy chain), length 115-140
            # VL: starts with DIQ/DIV/DTV/EIV (light chain kappa/lambda), length 100-120
            # Target: anything else (usually shorter, different sequence patterns)
            vh_sequence = ""
            vl_sequence = ""

            vh_patterns = ("EVQL", "QVQL", "DVQL", "EVKL", "QVKL", "EVLL")
            vl_patterns = ("DIQ", "DIV", "DTV", "EIV", "AIQ", "SIV", "DIQM", "DIQLT", "EIVL", "EIVMT")

            for entity_id, seq in sequences_by_entity.items():
                seq_len = len(seq)
                seq_prefix = seq[:5].upper()

                # Check VH patterns first (heavy chain variable region)
                if 115 <= seq_len <= 145 and any(seq.upper().startswith(p) for p in vh_patterns):
                    if not vh_sequence or seq_len > len(vh_sequence):
                        vh_sequence = seq
                        print(f"    -> Entity {entity_id} identified as VH ({seq_len} aa)")
                # Check VL patterns (light chain variable region)
                elif 95 <= seq_len <= 125 and any(seq.upper().startswith(p) for p in vl_patterns):
                    if not vl_sequence or seq_len > len(vl_sequence):
                        vl_sequence = seq
                        print(f"    -> Entity {entity_id} identified as VL ({seq_len} aa)")

            # If pattern matching didn't work, fall back to length heuristics
            if not vh_sequence or not vl_sequence:
                candidates = sorted(sequences_by_entity.items(), key=lambda x: len(x[1]), reverse=True)
                for entity_id, seq in candidates:
                    seq_len = len(seq)
                    if not vh_sequence and 115 <= seq_len <= 145:
                        vh_sequence = seq
                        print(f"    -> Entity {entity_id} assigned as VH by length ({seq_len} aa)")
                    elif not vl_sequence and 95 <= seq_len <= 125 and seq != vh_sequence:
                        vl_sequence = seq
                        print(f"    -> Entity {entity_id} assigned as VL by length ({seq_len} aa)")

            if vh_sequence and vl_sequence:
                design = {
                    "vh_sequence": vh_sequence,
                    "vl_sequence": vl_sequence,
                    "sequence": vh_sequence,  # For compatibility
                    "header": cif_file.stem,
                    "source_file": cif_file.name,
                    "binder_type": "fab",
                }
                designs.append(design)
                print(f"    -> Extracted VH ({len(vh_sequence)} aa) + VL ({len(vl_sequence)} aa)")
            else:
                print(f"    -> Could not identify VH/VL (VH={len(vh_sequence) if vh_sequence else 0}, VL={len(vl_sequence) if vl_sequence else 0})")

        except Exception as e:
            print(f"    Error parsing {cif_file.name}: {e}")
            import traceback
            traceback.print_exc()

    # Add metrics to designs
    for i, design in enumerate(designs):
        design["design_idx"] = i
        design_name = design.get("header", "")
        if design_name in metrics_by_design:
            metrics = metrics_by_design[design_name]
            try:
                design["ipTM"] = float(metrics.get("design_to_target_iptm", 0))
                design["pTM"] = float(metrics.get("design_ptm", 0))
                design["rmsd"] = float(metrics.get("filter_rmsd", 0))
                design["rmsd_design"] = float(metrics.get("filter_rmsd_design", 0))
                design["boltzgen_rank"] = metrics.get("_boltzgen_rank", i + 1)
            except (ValueError, TypeError):
                pass
        else:
            design["boltzgen_rank"] = i + 1

    print(f"Parsed {len(designs)} Fab designs")
    # Print RMSD summary
    rmsds = [d.get("rmsd", 0) for d in designs if d.get("rmsd")]
    if rmsds:
        print(f"RMSD range: {min(rmsds):.2f} - {max(rmsds):.2f} Å (threshold: 2.5 Å)")
    return designs


@app.function(
    image=boltzgen_image,
    volumes={models_dir: boltzgen_model_volume},
    timeout=120 * MINUTES,
    gpu="A100",
)
def run_boltzgen_batch(
    target_pdbs: list[dict],  # [{"name": "1XIW", "content": "...", "chain": "A"}]
    num_designs_per_target: int = 100,
    binder_length: int = 120,
    protocol: str = "nanobody-anything",
    seed: int = 42,
) -> dict[str, list[dict]]:
    """Run BoltzGen on multiple targets.

    Args:
        target_pdbs: List of dicts with "name", "content", and optional "chain" keys.
        num_designs_per_target: Designs per target.
        binder_length: Length of binder to design.
        protocol: BoltzGen protocol.
        seed: Base random seed.

    Returns:
        Dictionary mapping target names to design lists.
    """
    results = {}

    for i, target in enumerate(target_pdbs):
        target_name = target["name"]
        target_content = target["content"]
        target_chain = target.get("chain", "A")

        print(f"\nProcessing target {i+1}/{len(target_pdbs)}: {target_name}")

        try:
            designs = run_boltzgen(
                target_pdb_content=target_content,
                target_chain=target_chain,
                num_designs=num_designs_per_target,
                binder_length=binder_length,
                protocol=protocol,
                seed=seed + i * 1000,
            )
            results[target_name] = designs
        except Exception as e:
            print(f"Error designing for {target_name}: {e}")
            results[target_name] = [{"error": str(e)}]

    return results


@app.local_entrypoint()
def main(
    target_pdb: str = None,
    target_chain: str = "A",
    num_designs: int = 10,
    binder_length: int = 120,
    protocol: str = "nanobody-anything",
    design_type: str = "vhh",
    scaffold_dir: str = "data/fab_scaffolds",
    scaffolds: str = "adalimumab",
    seed: int = 42,
    download: bool = False,
):
    """Local entrypoint for testing.

    Args:
        target_pdb: Path to target PDB file.
        target_chain: Chain ID of target to extract (default: "A").
        num_designs: Number of designs.
        binder_length: Length of binder (120 for VHH).
        protocol: BoltzGen protocol.
        design_type: Type of design: "vhh" or "fab".
        scaffold_dir: Directory containing Fab scaffold files.
        scaffolds: Comma-separated list of scaffold names for Fab design.
        seed: Random seed.
        download: If True, download model weights and exit.
    """
    if download:
        print("Downloading BoltzGen model weights...")
        download_model.remote()
        print("Done!")
        return

    if target_pdb is None:
        print("Usage: modal run modal/boltzgen_app.py --target-pdb <path> [--target-chain A]")
        print("   or: modal run modal/boltzgen_app.py --download")
        print("\nOptions:")
        print("  --num-designs N      Number of designs to generate (default: 10)")
        print("  --design-type TYPE   Design type: vhh or fab (default: vhh)")
        print("  --binder-length N    Binder length in aa (default: 120, VHH only)")
        print("  --protocol NAME      Protocol: nanobody-anything, antibody-anything")
        print("  --scaffold-dir DIR   Directory with Fab scaffold files (fab only)")
        print("  --scaffolds NAMES    Comma-separated scaffold names (fab only)")
        print("\nExamples:")
        print("  # VHH design:")
        print("  modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb --design-type vhh")
        print("")
        print("  # Fab CDR redesign:")
        print("  modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb --design-type fab --scaffolds adalimumab,belimumab")
        return

    # Read target PDB
    with open(target_pdb, "r") as f:
        target_content = f.read()

    print(f"Running BoltzGen on {target_pdb}")
    print(f"  Target chain: {target_chain}")
    print(f"  Design type: {design_type}")
    print(f"  Num designs: {num_designs}")
    print(f"  Seed: {seed}")

    if design_type == "fab":
        # Fab CDR redesign
        scaffold_names = [s.strip() for s in scaffolds.split(",")]
        print(f"  Scaffolds: {scaffold_names}")

        # Load scaffold files
        scaffold_yaml_contents = {}
        scaffold_cif_contents = {}
        scaffold_path = Path(scaffold_dir)

        for name in scaffold_names:
            # Find YAML file
            yaml_files = list(scaffold_path.glob(f"{name}*.yaml"))
            if not yaml_files:
                print(f"ERROR: No YAML file found for scaffold '{name}' in {scaffold_dir}")
                return
            yaml_file = yaml_files[0]

            # Find CIF file
            cif_files = list(scaffold_path.glob(f"{name}*.cif"))
            if not cif_files:
                print(f"ERROR: No CIF file found for scaffold '{name}' in {scaffold_dir}")
                return
            cif_file = cif_files[0]

            scaffold_yaml_contents[name] = yaml_file.read_text()
            scaffold_cif_contents[name] = cif_file.read_text()
            print(f"  Loaded scaffold: {name} ({yaml_file.name}, {cif_file.name})")

        # Run Fab design on Modal
        designs = run_boltzgen_fab.remote(
            target_pdb_content=target_content,
            target_chain=target_chain,
            scaffold_names=scaffold_names,
            scaffold_yaml_contents=scaffold_yaml_contents,
            scaffold_cif_contents=scaffold_cif_contents,
            num_designs=num_designs,
            seed=seed,
        )

        print(f"\nGenerated {len(designs)} Fab designs:")
        for i, d in enumerate(designs[:5]):
            vh = d.get("vh_sequence", "N/A")
            vl = d.get("vl_sequence", "N/A")
            print(f"  {i+1}. VH: {vh[:30]}... ({len(vh)} aa)")
            print(f"      VL: {vl[:30]}... ({len(vl)} aa)")

    else:
        # VHH design (existing code)
        print(f"  Protocol: {protocol}")
        print(f"  Binder length: {binder_length} aa")

        designs = run_boltzgen.remote(
            target_pdb_content=target_content,
            target_chain=target_chain,
            num_designs=num_designs,
            binder_length=binder_length,
            protocol=protocol,
            seed=seed,
        )

        print(f"\nGenerated {len(designs)} designs:")
        for i, d in enumerate(designs[:5]):
            seq = d.get("sequence", "N/A")
            conf = d.get("confidence", "N/A")
            print(f"  {i+1}. Conf: {conf}, Seq: {seq[:40]}...")

    if len(designs) > 5:
        print(f"  ... and {len(designs) - 5} more")

    # Save results
    output_path = Path("boltzgen_designs.json")
    output_path.write_text(json.dumps(designs, indent=2))
    print(f"\nSaved to {output_path}")
