"""PDB file parsing and manipulation utilities."""

from pathlib import Path
from typing import Optional
import re


# Standard amino acid 3-letter to 1-letter mapping
AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def extract_sequence_from_pdb(
    pdb_path: str,
    chain_id: str = "A",
) -> str:
    """Extract amino acid sequence from a PDB file.

    Args:
        pdb_path: Path to PDB file.
        chain_id: Chain ID to extract.

    Returns:
        Amino acid sequence (1-letter codes).
    """
    sequence = []
    seen_residues = set()

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                if chain != chain_id:
                    continue

                res_name = line[17:20].strip()
                res_num = int(line[22:26])
                res_key = (res_num, res_name)

                if res_key not in seen_residues and res_name in AA_3TO1:
                    seen_residues.add(res_key)
                    sequence.append((res_num, AA_3TO1[res_name]))

    # Sort by residue number and extract sequence
    sequence.sort(key=lambda x: x[0])
    return "".join([aa for _, aa in sequence])


def extract_chain(
    pdb_path: str,
    chain_id: str,
    output_path: Optional[str] = None,
) -> str:
    """Extract a single chain from a PDB file.

    Args:
        pdb_path: Path to input PDB file.
        chain_id: Chain ID to extract.
        output_path: Path to save extracted chain (optional).

    Returns:
        PDB string for the extracted chain.
    """
    lines = []

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[21] == chain_id:
                    lines.append(line)
            elif line.startswith("TER"):
                if lines and lines[-1][21] == chain_id:
                    lines.append(line)
            elif line.startswith("END"):
                lines.append(line)
                break

    pdb_string = "".join(lines)

    if output_path:
        with open(output_path, "w") as f:
            f.write(pdb_string)

    return pdb_string


def get_chain_ids(pdb_path: str) -> list[str]:
    """Get all chain IDs in a PDB file.

    Args:
        pdb_path: Path to PDB file.

    Returns:
        List of chain IDs.
    """
    chain_ids = set()

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                chain_ids.add(line[21])

    return sorted(chain_ids)


def get_residue_positions(
    pdb_path: str,
    chain_id: str = "A",
) -> list[int]:
    """Get list of residue positions in a chain.

    Args:
        pdb_path: Path to PDB file.
        chain_id: Chain ID.

    Returns:
        List of residue numbers.
    """
    positions = set()

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == chain_id:
                res_num = int(line[22:26])
                positions.add(res_num)

    return sorted(positions)


def calculate_interface_residues(
    pdb_string: str,
    chain_a: str,
    chain_b: str,
    distance_cutoff: float = 5.0,
) -> tuple[list[int], list[int]]:
    """Calculate interface residues between two chains.

    Args:
        pdb_string: PDB file content as string.
        chain_a: First chain ID.
        chain_b: Second chain ID.
        distance_cutoff: Distance cutoff for contacts (Å).

    Returns:
        Tuple of (chain_a_residues, chain_b_residues) at interface.
    """
    # Parse atom coordinates
    atoms_a = []
    atoms_b = []

    for line in pdb_string.split("\n"):
        if not line.startswith("ATOM"):
            continue

        chain = line[21]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        res_num = int(line[22:26])

        if chain == chain_a:
            atoms_a.append((res_num, x, y, z))
        elif chain == chain_b:
            atoms_b.append((res_num, x, y, z))

    # Find interface residues
    interface_a = set()
    interface_b = set()

    cutoff_sq = distance_cutoff ** 2

    for res_a, xa, ya, za in atoms_a:
        for res_b, xb, yb, zb in atoms_b:
            dist_sq = (xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2
            if dist_sq <= cutoff_sq:
                interface_a.add(res_a)
                interface_b.add(res_b)

    return sorted(interface_a), sorted(interface_b)


def count_contacts(
    pdb_string: str,
    chain_a: str,
    chain_b: str,
    distance_cutoff: float = 5.0,
) -> int:
    """Count residue-residue contacts between two chains.

    A contact is counted when any atom of residue A is within the
    distance cutoff of any atom of residue B. Each residue pair
    is counted only once.

    Args:
        pdb_string: PDB file content as string.
        chain_a: First chain ID.
        chain_b: Second chain ID.
        distance_cutoff: Distance cutoff for contacts (Å).

    Returns:
        Number of residue-residue contacts.
    """
    # Parse atom coordinates with residue info
    atoms_a = []  # (res_num, x, y, z)
    atoms_b = []

    for line in pdb_string.split("\n"):
        if not line.startswith("ATOM"):
            continue

        chain = line[21]
        res_num = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        if chain == chain_a:
            atoms_a.append((res_num, x, y, z))
        elif chain == chain_b:
            atoms_b.append((res_num, x, y, z))

    # Count unique residue pairs in contact
    contact_pairs = set()
    cutoff_sq = distance_cutoff ** 2

    for res_a, xa, ya, za in atoms_a:
        for res_b, xb, yb, zb in atoms_b:
            dist_sq = (xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2
            if dist_sq <= cutoff_sq:
                contact_pairs.add((res_a, res_b))

    return len(contact_pairs)


def estimate_interface_area(
    pdb_string: str,
    chain_a: str,
    chain_b: str,
) -> float:
    """Estimate buried surface area at interface.

    This is a simplified estimation based on number of
    interface residues. For accurate BSA calculation,
    use specialized tools like FreeSASA.

    Args:
        pdb_string: PDB file content as string.
        chain_a: First chain ID.
        chain_b: Second chain ID.

    Returns:
        Estimated interface area in Å².
    """
    interface_a, interface_b = calculate_interface_residues(
        pdb_string, chain_a, chain_b
    )

    # Rough estimate: ~100 Å² per interface residue
    # This is a simplification - actual BSA varies by residue type
    total_interface_residues = len(interface_a) + len(interface_b)
    estimated_area = total_interface_residues * 50.0  # ~50 Å² per residue contribution

    return estimated_area


def renumber_chain(
    pdb_string: str,
    chain_id: str,
    start_number: int = 1,
) -> str:
    """Renumber residues in a chain starting from a given number.

    Args:
        pdb_string: PDB file content as string.
        chain_id: Chain to renumber.
        start_number: Starting residue number.

    Returns:
        PDB string with renumbered chain.
    """
    lines = []
    current_res = None
    new_res_num = start_number - 1

    for line in pdb_string.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if line[21] == chain_id:
                res_num = line[22:27].strip()
                if res_num != current_res:
                    current_res = res_num
                    new_res_num += 1

                # Replace residue number
                new_line = line[:22] + f"{new_res_num:4d} " + line[27:]
                lines.append(new_line)
            else:
                lines.append(line)
        else:
            lines.append(line)

    return "\n".join(lines)


def combine_pdbs(
    pdb_strings: list[str],
    chain_ids: Optional[list[str]] = None,
) -> str:
    """Combine multiple PDB strings into one.

    Args:
        pdb_strings: List of PDB content strings.
        chain_ids: Optional chain IDs to assign (default: A, B, C, ...).

    Returns:
        Combined PDB string.
    """
    if chain_ids is None:
        chain_ids = [chr(ord("A") + i) for i in range(len(pdb_strings))]

    combined_lines = []

    for pdb_string, chain_id in zip(pdb_strings, chain_ids):
        for line in pdb_string.split("\n"):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Update chain ID
                new_line = line[:21] + chain_id + line[22:]
                combined_lines.append(new_line)
            elif line.startswith("TER"):
                combined_lines.append(f"TER   {' ' * 73}\n")

    combined_lines.append("END\n")

    return "\n".join(combined_lines)


def extract_epitope_from_complex(
    pdb_string: str,
    target_chain: str,
    binder_chains: list[str],
    distance_cutoff: float = 5.0,
) -> list[int]:
    """Extract epitope residues on target chain that contact binder chains.

    This function identifies residues on the target chain that are within
    the distance cutoff of any atom in the binder chain(s).

    Args:
        pdb_string: PDB file content as string.
        target_chain: Chain ID of the target (e.g., "E" for CD3ε).
        binder_chains: Chain IDs of the binder (e.g., ["H", "L"] for Fab).
        distance_cutoff: Distance cutoff for contacts (Å).

    Returns:
        Sorted list of target residue numbers that form the epitope.
    """
    # Parse atom coordinates
    target_atoms = []  # (res_num, x, y, z)
    binder_atoms = []  # (x, y, z)

    for line in pdb_string.split("\n"):
        if not line.startswith("ATOM"):
            continue

        chain = line[21]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        res_num = int(line[22:26])

        if chain == target_chain:
            target_atoms.append((res_num, x, y, z))
        elif chain in binder_chains:
            binder_atoms.append((x, y, z))

    if not target_atoms or not binder_atoms:
        return []

    # Find target residues in contact with binder
    epitope_residues = set()
    cutoff_sq = distance_cutoff ** 2

    for res_num, xt, yt, zt in target_atoms:
        for xb, yb, zb in binder_atoms:
            dist_sq = (xt - xb) ** 2 + (yt - yb) ** 2 + (zt - zb) ** 2
            if dist_sq <= cutoff_sq:
                epitope_residues.add(res_num)
                break  # Found a contact, no need to check more binder atoms

    return sorted(epitope_residues)


def download_pdb(pdb_id: str, output_path: Optional[str] = None) -> str:
    """Download a PDB file from RCSB.

    Args:
        pdb_id: 4-letter PDB ID (e.g., "1SY6").
        output_path: Optional path to save the file.

    Returns:
        PDB file content as string.

    Raises:
        RuntimeError: If download fails.
    """
    import urllib.request
    import urllib.error

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"

    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            pdb_content = response.read().decode("utf-8")
    except urllib.error.URLError as e:
        raise RuntimeError(f"Failed to download PDB {pdb_id}: {e}")

    if output_path:
        with open(output_path, "w") as f:
            f.write(pdb_content)

    return pdb_content


def get_okt3_epitope_from_1sy6(
    pdb_path: Optional[str] = None,
    distance_cutoff: float = 5.0,
    cache_dir: str = "data/targets",
) -> list[int]:
    """Extract OKT3 epitope residues on CD3ε from 1SY6 structure.

    This function dynamically extracts the epitope residues from the
    1SY6 crystal structure rather than using hardcoded values. This
    ensures the epitope definition matches the actual structural contacts.

    Chain assignments in 1SY6 (from RCSB):
    - Chain A: CD3ε (epsilon) - the target antigen
    - Chain H: OKT3 VH (heavy chain variable region)
    - Chain L: OKT3 VL (light chain variable region)

    Note: Some documentation may refer to CD3ε as chain E, but the actual
    PDB file uses chain A.

    Args:
        pdb_path: Path to 1SY6 PDB file. If None, will try to load from
            cache_dir or download from RCSB.
        distance_cutoff: Distance cutoff for contacts (Å).
        cache_dir: Directory to cache downloaded PDB files.

    Returns:
        Sorted list of CD3ε residue numbers that form the OKT3 epitope.
    """
    import os

    pdb_content = None

    # Try to load from provided path
    if pdb_path and os.path.exists(pdb_path):
        with open(pdb_path, "r") as f:
            pdb_content = f.read()

    # Try to load from cache directory
    if pdb_content is None:
        cached_path = os.path.join(cache_dir, "1SY6.pdb")
        if os.path.exists(cached_path):
            with open(cached_path, "r") as f:
                pdb_content = f.read()

    # Download if not available locally
    if pdb_content is None:
        os.makedirs(cache_dir, exist_ok=True)
        cached_path = os.path.join(cache_dir, "1SY6.pdb")
        pdb_content = download_pdb("1SY6", cached_path)

    # Extract epitope: CD3ε residues (chain A) contacting OKT3 Fab (chains H, L)
    epitope = extract_epitope_from_complex(
        pdb_content,
        target_chain="A",  # CD3ε is chain A in 1SY6
        binder_chains=["H", "L"],
        distance_cutoff=distance_cutoff,
    )

    return epitope


def align_residue_numbers(
    query_sequence: str,
    reference_sequence: str,
    reference_residue_numbers: list[int],
) -> list[int]:
    """Map residue numbers from reference to query sequence via alignment.

    This is useful when comparing epitopes between structures with
    different numbering schemes. Uses simple pairwise alignment.

    Args:
        query_sequence: Sequence to map residues to.
        reference_sequence: Sequence with known residue numbers.
        reference_residue_numbers: Residue numbers in reference to map.

    Returns:
        Corresponding residue numbers in query sequence (1-indexed).
        Returns empty list if alignment fails.
    """
    # Simple alignment: find matching positions
    # For production use, consider using BioPython's pairwise2 or parasail

    # Build reference position to residue mapping
    ref_positions = {i: num for i, num in enumerate(reference_residue_numbers)}

    # Try to find query positions that match reference residues
    # This is a simplified approach - assumes sequences are similar length
    query_numbers = []

    if len(query_sequence) == len(reference_sequence):
        # Same length: direct mapping
        for ref_num in reference_residue_numbers:
            # Find which position in reference corresponds to this number
            # Assuming 1-indexed numbering starting at some offset
            query_numbers.append(ref_num)
    else:
        # Different lengths: need actual alignment
        # For now, return empty and warn
        import warnings
        warnings.warn(
            f"Sequence length mismatch ({len(query_sequence)} vs {len(reference_sequence)}). "
            "Residue number mapping requires proper sequence alignment."
        )
        return []

    return query_numbers
