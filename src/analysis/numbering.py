"""Antibody numbering and CDR identification using ANARCI.

Provides consistent numbering (IMGT, Chothia, Kabat) and CDR extraction
for antibody variable region sequences.
"""

from dataclasses import dataclass
from typing import Optional
import warnings


@dataclass
class NumberedResidue:
    """A single numbered residue."""

    position: str  # e.g., "27", "111A"
    amino_acid: str
    region: str  # "FR1", "CDR1", "FR2", etc.


@dataclass
class NumberedSequence:
    """A numbered antibody sequence with CDR annotations."""

    sequence: str
    chain_type: str  # "H" or "L"
    scheme: str  # "imgt", "chothia", "kabat"
    residues: list[NumberedResidue]
    cdr1: str
    cdr2: str
    cdr3: str
    cdr1_positions: tuple[int, int]  # (start, end) in original sequence
    cdr2_positions: tuple[int, int]
    cdr3_positions: tuple[int, int]
    framework1: str
    framework2: str
    framework3: str
    framework4: str
    species: Optional[str] = None
    v_gene: Optional[str] = None

    def get_cdr_positions_dict(self) -> dict[str, tuple[int, int]]:
        """Get CDR positions as a dictionary."""
        prefix = self.chain_type
        return {
            f"{prefix}1": self.cdr1_positions,
            f"{prefix}2": self.cdr2_positions,
            f"{prefix}3": self.cdr3_positions,
        }


def number_sequence(
    sequence: str,
    scheme: str = "imgt",
    chain_type: Optional[str] = None,
) -> Optional[NumberedSequence]:
    """Number an antibody sequence using ANARCI.

    Args:
        sequence: Amino acid sequence (VH or VL).
        scheme: Numbering scheme ('imgt', 'chothia', 'kabat', 'martin').
        chain_type: Force chain type ('H' or 'L'). If None, auto-detect.

    Returns:
        NumberedSequence with CDR annotations, or None if numbering fails.
    """
    try:
        from anarci import anarci, number
    except ImportError:
        warnings.warn(
            "ANARCI not installed. Install with: pip install anarci"
        )
        return None

    # Run ANARCI
    results = anarci([("query", sequence)], scheme=scheme, output=False)

    if not results or not results[0] or not results[0][0]:
        return None

    # ANARCI returns a tuple of 3 lists:
    # results[0] = list of numbering results (one per input sequence)
    # results[1] = list of chain info (one per input sequence)
    # results[2] = list of alignment details
    #
    # For each sequence, numbering result is: (numbering_list, start_idx, end_idx)
    # Chain info is a list of dicts with chain_type, species, etc.

    numbering_data = results[0][0][0]  # First sequence, first domain, numbering tuple
    chain_info = results[1][0] if len(results) > 1 and results[1] else None  # Chain info list

    # Handle both old and new ANARCI output formats
    if isinstance(numbering_data, tuple) and len(numbering_data) >= 1:
        # Format: (numbering_list, start, end)
        numbering = numbering_data[0]
    else:
        # Direct numbering list
        numbering = numbering_data

    # Extract chain type from chain_info
    # chain_info is a list of dicts, one per domain found
    detected_chain = "H"
    species = None
    v_gene = None
    if chain_info and len(chain_info) > 0 and isinstance(chain_info[0], dict):
        detected_chain = chain_info[0].get("chain_type", "H")
        species = chain_info[0].get("species")
        v_gene = chain_info[0].get("v_gene")

    chain = chain_type or detected_chain

    # Define CDR boundaries by scheme
    if scheme == "imgt":
        cdr_bounds = {
            "H": {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)},
            "L": {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)},
        }
    elif scheme == "chothia":
        cdr_bounds = {
            "H": {"CDR1": (26, 32), "CDR2": (52, 56), "CDR3": (95, 102)},
            "L": {"CDR1": (24, 34), "CDR2": (50, 52), "CDR3": (89, 97)},  # CDR-L2 is only 3 residues in Chothia
        }
    elif scheme == "kabat":
        cdr_bounds = {
            "H": {"CDR1": (31, 35), "CDR2": (50, 65), "CDR3": (95, 102)},
            "L": {"CDR1": (24, 34), "CDR2": (50, 56), "CDR3": (89, 97)},
        }
    else:
        cdr_bounds = {
            "H": {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)},
            "L": {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)},
        }

    bounds = cdr_bounds.get(chain, cdr_bounds["H"])

    # Extract regions
    residues = []
    cdr1_seq, cdr2_seq, cdr3_seq = [], [], []
    fr1_seq, fr2_seq, fr3_seq, fr4_seq = [], [], [], []
    cdr1_pos, cdr2_pos, cdr3_pos = [None, None], [None, None], [None, None]

    seq_idx = 0
    for pos, aa in numbering:
        if aa == "-":
            continue

        pos_num = pos[0]
        region = _get_region(pos_num, bounds)

        residues.append(NumberedResidue(
            position=f"{pos[0]}{pos[1] if pos[1] != ' ' else ''}",
            amino_acid=aa,
            region=region,
        ))

        # Collect sequences by region
        if region == "CDR1":
            cdr1_seq.append(aa)
            if cdr1_pos[0] is None:
                cdr1_pos[0] = seq_idx
            cdr1_pos[1] = seq_idx
        elif region == "CDR2":
            cdr2_seq.append(aa)
            if cdr2_pos[0] is None:
                cdr2_pos[0] = seq_idx
            cdr2_pos[1] = seq_idx
        elif region == "CDR3":
            cdr3_seq.append(aa)
            if cdr3_pos[0] is None:
                cdr3_pos[0] = seq_idx
            cdr3_pos[1] = seq_idx
        elif region == "FR1":
            fr1_seq.append(aa)
        elif region == "FR2":
            fr2_seq.append(aa)
        elif region == "FR3":
            fr3_seq.append(aa)
        elif region == "FR4":
            fr4_seq.append(aa)

        seq_idx += 1

    return NumberedSequence(
        sequence=sequence,
        chain_type=chain,
        scheme=scheme,
        residues=residues,
        cdr1="".join(cdr1_seq),
        cdr2="".join(cdr2_seq),
        cdr3="".join(cdr3_seq),
        cdr1_positions=tuple(cdr1_pos) if cdr1_pos[0] is not None else (0, 0),
        cdr2_positions=tuple(cdr2_pos) if cdr2_pos[0] is not None else (0, 0),
        cdr3_positions=tuple(cdr3_pos) if cdr3_pos[0] is not None else (0, 0),
        framework1="".join(fr1_seq),
        framework2="".join(fr2_seq),
        framework3="".join(fr3_seq),
        framework4="".join(fr4_seq),
        species=species,
        v_gene=v_gene,
    )


def _get_region(pos_num: int, bounds: dict) -> str:
    """Determine which region a position falls into."""
    cdr1 = bounds["CDR1"]
    cdr2 = bounds["CDR2"]
    cdr3 = bounds["CDR3"]

    if pos_num < cdr1[0]:
        return "FR1"
    elif cdr1[0] <= pos_num <= cdr1[1]:
        return "CDR1"
    elif cdr1[1] < pos_num < cdr2[0]:
        return "FR2"
    elif cdr2[0] <= pos_num <= cdr2[1]:
        return "CDR2"
    elif cdr2[1] < pos_num < cdr3[0]:
        return "FR3"
    elif cdr3[0] <= pos_num <= cdr3[1]:
        return "CDR3"
    else:
        return "FR4"


def get_cdr_positions(
    sequence: str,
    chain_type: str = "H",
    scheme: str = "imgt",
) -> dict[str, tuple[int, int]]:
    """Get CDR positions for a sequence.

    Args:
        sequence: Amino acid sequence.
        chain_type: 'H' for heavy, 'L' for light.
        scheme: Numbering scheme.

    Returns:
        Dict mapping CDR names to (start, end) positions in the sequence.
    """
    numbered = number_sequence(sequence, scheme=scheme, chain_type=chain_type)

    if numbered is None:
        # Return empty positions if numbering fails
        prefix = chain_type
        return {
            f"{prefix}1": (0, 0),
            f"{prefix}2": (0, 0),
            f"{prefix}3": (0, 0),
        }

    return numbered.get_cdr_positions_dict()


def extract_cdrs(
    sequence: str,
    chain_type: str = "H",
    scheme: str = "imgt",
) -> dict[str, str]:
    """Extract CDR sequences from an antibody variable region.

    Args:
        sequence: Amino acid sequence (VH or VL).
        chain_type: 'H' for heavy, 'L' for light.
        scheme: Numbering scheme.

    Returns:
        Dict with CDR1, CDR2, CDR3 sequences.
    """
    numbered = number_sequence(sequence, scheme=scheme, chain_type=chain_type)

    if numbered is None:
        return {"CDR1": "", "CDR2": "", "CDR3": ""}

    return {
        "CDR1": numbered.cdr1,
        "CDR2": numbered.cdr2,
        "CDR3": numbered.cdr3,
    }


def extract_frameworks(
    sequence: str,
    chain_type: str = "H",
    scheme: str = "imgt",
) -> dict[str, str]:
    """Extract framework sequences from an antibody variable region.

    Args:
        sequence: Amino acid sequence (VH or VL).
        chain_type: 'H' for heavy, 'L' for light.
        scheme: Numbering scheme.

    Returns:
        Dict with FR1, FR2, FR3, FR4 sequences.
    """
    numbered = number_sequence(sequence, scheme=scheme, chain_type=chain_type)

    if numbered is None:
        return {"FR1": "", "FR2": "", "FR3": "", "FR4": ""}

    return {
        "FR1": numbered.framework1,
        "FR2": numbered.framework2,
        "FR3": numbered.framework3,
        "FR4": numbered.framework4,
    }


def is_valid_antibody_sequence(
    sequence: str,
    chain_type: Optional[str] = None,
) -> bool:
    """Check if a sequence can be numbered as an antibody.

    Args:
        sequence: Amino acid sequence to check.
        chain_type: Expected chain type, or None to accept any.

    Returns:
        True if sequence is a valid antibody variable region.
    """
    numbered = number_sequence(sequence)

    if numbered is None:
        return False

    if chain_type and numbered.chain_type != chain_type:
        return False

    # Check that CDRs were found
    if not numbered.cdr1 or not numbered.cdr3:
        return False

    return True


def get_imgt_to_sequence_mapping(
    sequence: str,
    chain_type: str = "H",
) -> dict[int, int]:
    """Create a mapping from IMGT position numbers to sequence indices.

    This is critical for applying mutations specified in IMGT numbering
    to raw sequences. IMGT numbering has gaps and insertions that don't
    correspond directly to sequence indices.

    Args:
        sequence: Amino acid sequence (VH or VL).
        chain_type: 'H' for heavy, 'L' for light.

    Returns:
        Dict mapping IMGT position (int) to 0-indexed sequence position.
        For positions with insertion codes (e.g., 111A), only the base
        position is included in the mapping.

    Example:
        >>> mapping = get_imgt_to_sequence_mapping(vh_sequence, "H")
        >>> seq_idx = mapping.get(50)  # Get sequence index for IMGT position 50
        >>> if seq_idx is not None:
        ...     residue = vh_sequence[seq_idx]
    """
    numbered = number_sequence(sequence, scheme="imgt", chain_type=chain_type)

    if numbered is None:
        return {}

    mapping = {}
    for idx, residue in enumerate(numbered.residues):
        # Parse position string (may have insertion code like "111A")
        pos_str = residue.position
        # Extract base position number (ignore insertion codes for mapping)
        base_pos = int("".join(c for c in pos_str if c.isdigit()))

        # Only store first occurrence (base position without insertion)
        if base_pos not in mapping:
            mapping[base_pos] = idx

    return mapping


def get_sequence_to_imgt_mapping(
    sequence: str,
    chain_type: str = "H",
) -> dict[int, str]:
    """Create a mapping from sequence indices to IMGT position strings.

    Args:
        sequence: Amino acid sequence (VH or VL).
        chain_type: 'H' for heavy, 'L' for light.

    Returns:
        Dict mapping 0-indexed sequence position to IMGT position string.
        Position strings may include insertion codes (e.g., "111", "111A").
    """
    numbered = number_sequence(sequence, scheme="imgt", chain_type=chain_type)

    if numbered is None:
        return {}

    return {idx: residue.position for idx, residue in enumerate(numbered.residues)}
