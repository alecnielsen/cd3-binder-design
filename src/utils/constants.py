"""Constants for antibody analysis and engineering."""

# Standard amino acids
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Amino acid properties
HYDROPHOBIC = set("AILMFVPWG")
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
DEAMIDATION_MOTIFS = ["NG", "NS", "NT", "ND"]  # N followed by small residue
ISOMERIZATION_MOTIFS = ["DG", "DS", "DT", "DD"]  # D followed by small residue
OXIDATION_RESIDUES = ["M", "W"]  # Methionine, Tryptophan susceptible

# Glycosylation pattern: N-X-S/T where X != P
GLYCOSYLATION_PATTERN = r"N[^P][ST]"

# CDR approximate positions (Chothia numbering)
CDR_RANGES_CHOTHIA = {
    "H1": (26, 32),
    "H2": (52, 56),
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

# Common scFv linker patterns for parsing
SCFV_LINKER_PATTERNS = [
    "GGGGSGGGGSGGGGS",  # (G4S)₃ - most common
    "GGGGSGGGGSGGGGSGGGGS",  # (G4S)₄
    "GGGGSGGGGS",  # (G4S)₂
    "GSTSGSGKPGSGEGSTKG",  # Whitlow linker
]


def parse_scfv(sequence: str, linker: str = None) -> tuple[str, str] | None:
    """Parse an scFv sequence into VH and VL components.

    Args:
        sequence: Full scFv sequence (VH-linker-VL format).
        linker: Specific linker to search for. If None, tries common linkers.

    Returns:
        Tuple of (vh, vl) if linker found, None otherwise.
    """
    if linker:
        linkers_to_try = [linker]
    else:
        linkers_to_try = SCFV_LINKER_PATTERNS

    for lnk in linkers_to_try:
        if lnk in sequence:
            idx = sequence.find(lnk)
            vh = sequence[:idx]
            vl = sequence[idx + len(lnk):]
            # Sanity check: VH and VL should be reasonable lengths
            if 100 <= len(vh) <= 140 and 100 <= len(vl) <= 130:
                return vh, vl

    return None


def is_likely_scfv(sequence: str) -> bool:
    """Check if a sequence is likely an scFv (VH-linker-VL).

    Args:
        sequence: Amino acid sequence to check.

    Returns:
        True if sequence appears to be scFv format.
    """
    # scFv typically 240-280 aa (VH ~120 + linker ~15 + VL ~110)
    if not (230 <= len(sequence) <= 300):
        return False

    # Check for common linker patterns
    for linker in SCFV_LINKER_PATTERNS:
        if linker in sequence:
            return True

    return False
