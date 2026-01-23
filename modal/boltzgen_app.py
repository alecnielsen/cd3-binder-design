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

import json
import subprocess
from pathlib import Path

import modal

MINUTES = 60

app = modal.App("boltzgen-cd3")

# Container with BoltzGen
boltzgen_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("boltzgen")
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
    target_sequence: str,
    binder_length: int = 120,
    hotspot_residues: list[int] | None = None,
) -> str:
    """Build BoltzGen design specification YAML.

    Args:
        target_sequence: Target protein sequence.
        binder_length: Length of binder to design (120 for VHH).
        hotspot_residues: Optional list of target residues to focus binding on.

    Returns:
        YAML string for BoltzGen design spec.
    """
    # Build design spec
    # Target is chain A (fixed), binder is chain B (to be designed)
    spec = {
        "entities": [
            {
                "protein": {
                    "id": "A",
                    "sequence": target_sequence,
                    "design": False,
                }
            },
            {
                "protein": {
                    "id": "B",
                    "sequence": "X" * binder_length,  # X marks designable positions
                    "design": True,
                    "binding": ["A"],
                }
            }
        ]
    }

    # Add hotspot constraints if specified
    if hotspot_residues:
        spec["entities"][1]["protein"]["constraints"] = {
            "contacts": [{"chain": "A", "residues": hotspot_residues}]
        }

    # Convert to YAML manually (avoid PyYAML dependency)
    lines = ["entities:"]

    for entity in spec["entities"]:
        lines.append("  - protein:")
        protein = entity["protein"]
        lines.append(f"      id: {protein['id']}")
        lines.append(f"      sequence: {protein['sequence']}")
        if "design" in protein:
            lines.append(f"      design: {'true' if protein['design'] else 'false'}")
        if "binding" in protein:
            lines.append(f"      binding: [{', '.join(protein['binding'])}]")
        if "constraints" in protein:
            lines.append("      constraints:")
            lines.append("        contacts:")
            for contact in protein["constraints"]["contacts"]:
                lines.append(f"          - chain: {contact['chain']}")
                residues_str = ", ".join(str(r) for r in contact["residues"])
                lines.append(f"            residues: [{residues_str}]")

    return "\n".join(lines)


@app.function(
    volumes={models_dir: boltzgen_model_volume},
    timeout=20 * MINUTES,
    image=download_image,
)
def download_model(force_download: bool = False):
    """Download BoltzGen model weights to Modal volume."""
    from huggingface_hub import snapshot_download

    # BoltzGen models - adjust repo_id if different
    snapshot_download(
        repo_id="boltz-community/boltzgen",
        local_dir=models_dir,
        force_download=force_download,
    )
    boltzgen_model_volume.commit()
    print(f"Model downloaded to {models_dir}")


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

    # Get target sequence
    target_sequence = parse_pdb_sequence(extracted_pdb, target_chain)
    print(f"Target sequence ({len(target_sequence)} aa): {target_sequence[:50]}...")

    # Build design spec YAML
    spec_yaml = build_design_spec_yaml(
        target_sequence=target_sequence,
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

    if result.returncode != 0:
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        raise RuntimeError(f"BoltzGen failed: {result.stderr}")

    print(f"STDOUT: {result.stdout}")

    # Parse outputs - BoltzGen outputs FASTA and/or CIF files
    designs = []

    # Look for output files
    fasta_files = list(output_dir.glob("**/*.fasta")) + list(output_dir.glob("**/*.fa"))
    cif_files = list(output_dir.glob("**/*.cif"))
    json_files = list(output_dir.glob("**/*.json"))

    print(f"Found {len(fasta_files)} FASTA, {len(cif_files)} CIF, {len(json_files)} JSON files")

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

    # Add design indices
    for i, design in enumerate(designs):
        design["design_idx"] = i

    print(f"Parsed {len(designs)} designs")
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
        print("  --binder-length N    Binder length in aa (default: 120 for VHH)")
        print("  --protocol NAME      Protocol: nanobody-anything, protein-anything")
        return

    # Read target PDB
    with open(target_pdb, "r") as f:
        target_content = f.read()

    print(f"Running BoltzGen on {target_pdb}")
    print(f"  Target chain: {target_chain}")
    print(f"  Protocol: {protocol}")
    print(f"  Binder length: {binder_length} aa")
    print(f"  Num designs: {num_designs}")
    print(f"  Seed: {seed}")

    # Run on Modal
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
