"""Constants for antibody analysis and engineering."""

from typing import Optional, Tuple

# Standard amino acids
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Amino acid properties
# Note: HYDROPHOBIC excludes G (minimal hydrophobicity, Kyte-Doolittle -0.4)
# and P (polar, Kyte-Doolittle -1.6) which caused false-positive patch detection
HYDROPHOBIC = set("AILMFVW")  # Aliphatic + aromatic (excluding polar G, P)
POLAR = set("STNQ")
CHARGED_POS = set("RKH")
CHARGED_NEG = set("DE")
AROMATIC = set("FWY")
SMALL = set("AGST")

# Kyte-Doolittle hydrophobicity scale
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Sequence liability motifs
# Deamidation: N followed by small/flexible residue enables backbone rotation
# NG is fastest, NS medium-high, NT lower, ND/NH moderate risk
# Reference: Robinson & Robinson (2001) PNAS 98:944-949
DEAMIDATION_MOTIFS = ["NG", "NS", "NT", "ND", "NH"]

# Isomerization: D followed by small/flexible residue enables succinimide formation
# Reference: Cacia et al. (1996) Biochemistry 35:1897-1903
ISOMERIZATION_MOTIFS = ["DG", "DS", "DT", "DD", "DH", "DN"]

# Oxidation: Methionine and Tryptophan are most susceptible
# Note: Histidine can also oxidize but at lower rate
OXIDATION_RESIDUES = ["M", "W"]  # Methionine, Tryptophan susceptible

# Glycosylation pattern: N-X-S/T where X != P
GLYCOSYLATION_PATTERN = r"N[^P][ST]"

# CDR approximate positions (Extended Chothia/AbM numbering)
# NOTE: Original Chothia uses H2: 52-56, L2: 50-52. We use extended definitions
# (H2: 52-58, L2: 50-56) which are more common in antibody engineering tools.
CDR_RANGES_CHOTHIA = {
    "H1": (26, 32),
    "H2": (52, 58),
    "H3": (95, 102),
    "L1": (24, 34),
    "L2": (50, 56),
    "L3": (89, 97),
}

# CDR positions (IMGT numbering)
CDR_RANGES_IMGT = {
    "H1": (27, 38),
    "H2": (56, 65),
    "H3": (105, 117),
    "L1": (27, 38),
    "L2": (56, 65),
    "L3": (105, 117),
}

# Linker sequences
LINKERS = {
    "G4S": "GGGGS",
    "G4S_2": "GGGGSGGGGS",
    "G4S_3": "GGGGSGGGGSGGGGS",
    "G4S_4": "GGGGSGGGGSGGGGSGGGGS",
    "WHITLOW": "GSTSGSGKPGSGEGSTKG",
}

# Knob-in-hole mutations (EU numbering)
KNOB_MUTATIONS = {"T366W": "W"}
HOLE_MUTATIONS = {"T366S": "S", "L368A": "A", "Y407V": "V"}

# LALA mutations for Fc silencing (EU numbering)
LALA_MUTATIONS = {"L234A": "A", "L235A": "A"}

# Default filtering thresholds
DEFAULT_FILTER_THRESHOLDS = {
    "min_pdockq": 0.5,
    "min_interface_area": 800,  # Angstrom^2
    "min_contacts": 10,
    "min_oasis_score": 0.8,
    "cdr_h3_length_min": 8,
    "cdr_h3_length_max": 20,
    "net_charge_min": -2,
    "net_charge_max": 4,
    "pi_min": 6.0,
    "pi_max": 9.0,
    "max_hydrophobic_patches": 2,
    "max_cdr_aromatic_fraction": 0.2,
}

# Common scFv linker patterns for parsing (ordered by specificity/length)
SCFV_LINKER_PATTERNS = [
    # Long G4S linkers (most common, check first)
    "GGGGSGGGGSGGGGSGGGGS",  # (G4S)₄
    "GGGGSGGGGSGGGGS",  # (G4S)₃ - most common
    "GGGGSGGGGS",  # (G4S)₂
    "GGGGS",  # (G4S)₁
    # Whitlow linker (218 design)
    "GSTSGSGKPGSGEGSTKG",
    # G3S linkers
    "GGGSGGGSGGGSGGGS",  # (G3S)₄
    "GGGSGGGSGGGSGGS",  # (G3S)₄ variant
    "GGGSGGGSGGS",  # (G3S)₃
    "GGGSGGGSGGGSGGGSGGS",  # (G3S)₅
    # Helical linkers
    "AEAAAKEAAAKEAAAKA",  # A(EAAAK)₃A
    "AEAAAKEAAAKA",  # A(EAAAK)₂A
    # Other common linkers
    "GSGSGSGSGS",  # (GS)₅
    "GGGSGGGS",  # (G3S)₂
    "KESGSVSSEQLAQFRSLD",  # Bird linker
    "EGKSSGSGSESKST",  # Alternative linker
]


def _find_linker_by_heuristic(sequence: str) -> Optional[Tuple[int, int]]:
    """Find a potential linker region using heuristics.

    Looks for glycine/serine-rich regions that could be linkers.

    Args:
        sequence: Full sequence to search.

    Returns:
        Tuple of (start_index, end_index) if found, None otherwise.
    """
    import re

    # Look for stretches of primarily G and S (at least 10 aa, >70% G/S)
    min_linker_len = 10
    max_linker_len = 30

    # Find all positions with G or S
    gs_positions = [i for i, aa in enumerate(sequence) if aa in "GS"]

    if len(gs_positions) < min_linker_len:
        return None

    # Sliding window to find dense G/S regions
    best_region = None
    best_density = 0

    for linker_len in range(min_linker_len, min(max_linker_len + 1, len(sequence) - 200)):
        for start in range(100, len(sequence) - 100 - linker_len):  # VH ~100-140, VL ~100-130
            end = start + linker_len
            region = sequence[start:end]
            gs_count = sum(1 for aa in region if aa in "GS")
            density = gs_count / linker_len

            # Must be >70% G/S and result in reasonable VH/VL lengths
            vh_len = start
            vl_len = len(sequence) - end

            if density > 0.7 and 100 <= vh_len <= 140 and 100 <= vl_len <= 130:
                # Prefer longer linkers and higher density
                score = density * linker_len
                if score > best_density:
                    best_density = score
                    best_region = (start, end)

    return best_region


def _score_scfv_split(vh_len: int, vl_len: int, linker_len: int) -> float:
    """Score a potential scFv split based on how typical the lengths are.

    Args:
        vh_len: Length of VH region.
        vl_len: Length of VL region.
        linker_len: Length of linker used.

    Returns:
        Score where higher is better. Returns -inf for invalid splits.
    """
    # Reject invalid lengths
    if not (100 <= vh_len <= 140 and 100 <= vl_len <= 130):
        return float("-inf")

    # Ideal lengths: VH ~120, VL ~110, linker ~15
    vh_score = -abs(vh_len - 120)  # Penalize deviation from 120
    vl_score = -abs(vl_len - 110)  # Penalize deviation from 110
    linker_score = linker_len * 0.5  # Prefer longer linkers (more specific)

    return vh_score + vl_score + linker_score


def parse_scfv(sequence: str, linker: str = None) -> Optional[Tuple[str, str]]:
    """Parse an scFv sequence into VH and VL components.

    Finds all occurrences of linker patterns and returns the split that
    produces the most plausible VH/VL lengths.

    Args:
        sequence: Full scFv sequence (VH-linker-VL format).
        linker: Specific linker to search for. If None, tries common linkers
            and falls back to heuristic detection.

    Returns:
        Tuple of (vh, vl) if linker found, None otherwise.
    """
    if linker:
        linkers_to_try = [linker]
    else:
        linkers_to_try = SCFV_LINKER_PATTERNS

    # Collect all valid splits with their scores
    candidates = []

    for lnk in linkers_to_try:
        # Find ALL occurrences of this linker pattern
        start = 0
        while True:
            idx = sequence.find(lnk, start)
            if idx == -1:
                break

            vh = sequence[:idx]
            vl = sequence[idx + len(lnk):]
            score = _score_scfv_split(len(vh), len(vl), len(lnk))

            if score > float("-inf"):
                candidates.append((vh, vl, score, lnk))

            start = idx + 1  # Continue searching for more occurrences

    # Return the best split if any valid candidates found
    if candidates:
        candidates.sort(key=lambda x: x[2], reverse=True)
        return candidates[0][0], candidates[0][1]

    # Fallback to heuristic detection if no explicit linker provided
    if linker is None:
        region = _find_linker_by_heuristic(sequence)
        if region:
            start, end = region
            vh = sequence[:start]
            vl = sequence[end:]
            if 100 <= len(vh) <= 140 and 100 <= len(vl) <= 130:
                return vh, vl

    return None


def is_likely_scfv(sequence: str) -> bool:
    """Check if a sequence is likely an scFv (VH-linker-VL).

    Uses both explicit linker pattern matching and heuristic detection.

    Args:
        sequence: Amino acid sequence to check.

    Returns:
        True if sequence appears to be scFv format.
    """
    # scFv typically 230-300 aa (VH ~120 + linker ~15-25 + VL ~110)
    if not (230 <= len(sequence) <= 320):
        return False

    # Check for explicit linker patterns
    for linker in SCFV_LINKER_PATTERNS:
        if linker in sequence:
            # Verify the split produces reasonable lengths
            idx = sequence.find(linker)
            vh_len = idx
            vl_len = len(sequence) - idx - len(linker)
            if 100 <= vh_len <= 140 and 100 <= vl_len <= 130:
                return True

    # Fallback: check if heuristic finds a plausible linker region
    region = _find_linker_by_heuristic(sequence)
    if region:
        return True

    return False
