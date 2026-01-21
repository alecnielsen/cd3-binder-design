"""BoltzGen Modal deployment for de novo binder design.

This module provides a Modal app for running BoltzGen
on GPU infrastructure for de novo antibody design.

IMPORTANT: This module extracts ONLY the target chain from input PDBs
before passing to BoltzGen. Multi-chain PDBs (e.g., with bound antibodies)
would bias design if passed in full.

Deploy with:
    modal deploy modal/boltzgen_app.py

Run locally:
    modal run modal/boltzgen_app.py::run_boltzgen --target-pdb data/targets/cd3.pdb
"""

import modal


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
        # List available chains for debugging
        available_chains = set()
        for line in pdb_content.split("\n"):
            if line.startswith("ATOM") and len(line) > 21:
                available_chains.add(line[21])
        raise ValueError(
            f"Chain '{chain_id}' not found in PDB. "
            f"Available chains: {sorted(available_chains)}"
        )

    return "\n".join(lines)


def validate_pdb_for_design(pdb_content: str, target_chain: str) -> tuple[bool, str]:
    """Validate that a PDB is appropriate for binder design.

    Checks for common issues that could bias design results:
    - Multi-chain complexes (e.g., target + bound antibody)
    - Known problematic structures (1XIW with UCHT1, 1SY6 with OKT3)

    Args:
        pdb_content: PDB file content as string.
        target_chain: Chain ID that will be used for design.

    Returns:
        Tuple of (is_valid, warning_message).
        is_valid is True if the PDB appears suitable for design.
        warning_message describes any issues found.
    """
    # Count chains in the PDB
    chains = set()
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") and len(line) > 21:
            chains.add(line[21])

    warnings = []

    # Check for multi-chain complexes
    if len(chains) > 1:
        other_chains = sorted(c for c in chains if c != target_chain)
        warnings.append(
            f"PDB contains {len(chains)} chains ({sorted(chains)}). "
            f"Only chain {target_chain} will be extracted for design. "
            f"Chains {other_chains} will be removed to prevent design bias."
        )

    # Check for known problematic structures
    # Look for UCHT1 or OKT3 signatures in HEADER/TITLE
    header_lines = [
        line for line in pdb_content.split("\n")[:50]
        if line.startswith(("HEADER", "TITLE", "COMPND"))
    ]
    header_text = " ".join(header_lines).upper()

    if "UCHT1" in header_text or "1XIW" in header_text:
        if len(chains) > 2:  # CD3ε + CD3δ + UCHT1
            warnings.append(
                "Detected 1XIW structure with UCHT1 Fab. "
                "Extracting only the target chain to avoid design bias."
            )

    if "OKT3" in header_text or "1SY6" in header_text:
        warnings.append(
            "Detected 1SY6 structure with OKT3 Fab. "
            "Note: Chain A is a CD3γ/ε fusion construct, not pure CD3ε. "
            "Extracting only the target chain to avoid design bias."
        )

    is_valid = True  # We always extract the single chain, so it's always "valid"
    warning_message = " ".join(warnings) if warnings else ""

    return is_valid, warning_message


# Create Modal app
app = modal.App("boltzgen-cd3")

# Define container image with BoltzGen dependencies
boltzgen_image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "torch>=2.0.0",
        "numpy",
        "biopython",
        "einops",
        "scipy",
    )
    # Note: BoltzGen installation may require additional steps
    # depending on the actual package distribution
    .pip_install("boltzgen")
)


@app.function(
    image=boltzgen_image,
    gpu="A100",  # Request A100 GPU
    timeout=3600,  # 1 hour timeout
    retries=2,
)
def run_boltzgen(
    target_pdb_content: str,
    target_chain: str = "A",
    binder_type: str = "vhh",
    num_designs: int = 100,
    hotspot_residues: list[int] = None,
    seed: int = 42,
    temperature: float = 1.0,
    num_recycles: int = 3,
) -> list[dict]:
    """Run BoltzGen to design binders.

    IMPORTANT: This function extracts ONLY the target chain before passing
    to BoltzGen. Multi-chain PDBs (e.g., CD3 + bound Fab) would bias designs
    if passed in full, potentially targeting the wrong epitope or generating
    antibody-like sequences instead of true CD3 binders.

    Args:
        target_pdb_content: PDB file content as string (can be multi-chain).
        target_chain: Chain ID of target in PDB (will be extracted).
        binder_type: Type of binder to design ("vhh" or "scfv").
        num_designs: Number of designs to generate.
        hotspot_residues: Optional residues to target for binding.
        seed: Random seed for reproducibility.
        temperature: Sampling temperature.
        num_recycles: Number of structure prediction recycles.

    Returns:
        List of design dictionaries with sequences and scores.

    Raises:
        ValueError: If target_chain is not found in the PDB.
    """
    import tempfile
    import os

    # CRITICAL: Extract only the target chain to avoid design bias
    # Multi-chain PDBs (e.g., 1XIW with UCHT1, 1SY6 with OKT3) would cause
    # BoltzGen to see the bound antibody, biasing designs inappropriately.
    extracted_pdb = extract_chain_from_pdb_content(target_pdb_content, target_chain)

    # Write extracted single-chain PDB to temp file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
        f.write(extracted_pdb)
        target_path = f.name

    try:
        # Import BoltzGen
        # Note: Actual import and API may differ based on BoltzGen version
        from boltzgen import BoltzGen

        # Load model
        model = BoltzGen.load()

        # Configure design parameters
        design_config = {
            "binder_type": binder_type,
            "num_samples": num_designs,
            "seed": seed,
            "temperature": temperature,
            "num_recycles": num_recycles,
        }

        if hotspot_residues:
            design_config["hotspot_residues"] = hotspot_residues

        # Run design
        designs = model.design(
            target_structure=target_path,
            target_chain=target_chain,
            **design_config,
        )

        # Convert to output format
        results = []
        for i, design in enumerate(designs):
            results.append({
                "sequence": design.sequence,
                "confidence": getattr(design, "confidence", 0.0),
                "plddt": getattr(design, "plddt", None),
                "ptm": getattr(design, "ptm", None),
                "design_idx": i,
            })

        return results

    finally:
        # Clean up temp file
        os.unlink(target_path)


@app.function(
    image=boltzgen_image,
    gpu="A100",
    timeout=7200,  # 2 hours for batch
)
def run_boltzgen_batch(
    target_pdbs: list[dict],  # [{"name": "1XIW", "content": "...", "chain": "A"}]
    binder_type: str = "vhh",
    num_designs_per_target: int = 100,
    seed: int = 42,
) -> dict[str, list[dict]]:
    """Run BoltzGen on multiple targets.

    IMPORTANT: Each target dict should include a "chain" key specifying which
    chain to extract for design. Only the target chain is passed to BoltzGen
    to avoid design bias from bound antibodies in crystal structures.

    Args:
        target_pdbs: List of dicts with "name", "content", and optional "chain" keys.
            If "chain" is not specified, defaults to "A".
        binder_type: Type of binder to design.
        num_designs_per_target: Designs per target.
        seed: Base random seed.

    Returns:
        Dictionary mapping target names to design lists.
    """
    results = {}

    for i, target in enumerate(target_pdbs):
        target_name = target["name"]
        target_content = target["content"]
        target_chain = target.get("chain", "A")  # Default to chain A

        designs = run_boltzgen.local(
            target_pdb_content=target_content,
            target_chain=target_chain,
            binder_type=binder_type,
            num_designs=num_designs_per_target,
            seed=seed + i * 1000,
        )

        results[target_name] = designs

    return results


@app.local_entrypoint()
def main(
    target_pdb: str = None,
    target_chain: str = "A",
    binder_type: str = "vhh",
    num_designs: int = 10,
    seed: int = 42,
):
    """Local entrypoint for testing.

    Args:
        target_pdb: Path to target PDB file.
        target_chain: Chain ID of target to extract (default: "A").
        binder_type: Type of binder ("vhh" or "scfv").
        num_designs: Number of designs.
        seed: Random seed.
    """
    if target_pdb is None:
        print("Usage: modal run modal/boltzgen_app.py --target-pdb <path> [--target-chain A]")
        return

    # Read target PDB
    with open(target_pdb, "r") as f:
        target_content = f.read()

    # Validate PDB and warn about potential issues
    is_valid, warning = validate_pdb_for_design(target_content, target_chain)
    if warning:
        print(f"\n⚠️  WARNING: {warning}\n")

    print(f"Running BoltzGen on {target_pdb}")
    print(f"  Target chain: {target_chain}")
    print(f"  Binder type: {binder_type}")
    print(f"  Num designs: {num_designs}")
    print(f"  Seed: {seed}")

    # Run on Modal (chain extraction happens inside run_boltzgen)
    designs = run_boltzgen.remote(
        target_pdb_content=target_content,
        target_chain=target_chain,
        binder_type=binder_type,
        num_designs=num_designs,
        seed=seed,
    )

    print(f"\nGenerated {len(designs)} designs:")
    for i, d in enumerate(designs[:5]):
        print(f"  {i+1}. Confidence: {d['confidence']:.3f}, Seq: {d['sequence'][:30]}...")

    if len(designs) > 5:
        print(f"  ... and {len(designs) - 5} more")
