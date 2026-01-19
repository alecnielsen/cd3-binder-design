"""Bispecific antibody formatting modules.

This package provides formatters for assembling CD3 binders into
various bispecific antibody architectures.

Available formats:
- CrossMab: Fab x Fab with CH1-CL swap for correct LC pairing
- Fab + scFv: Asymmetric with knob-in-hole
- Fab + VHH: Asymmetric with knob-in-hole (smallest asymmetric)
- IgG-scFv: Symmetric Morrison format (2x scFv on C-terminus)
- IgG-VHH: Symmetric Morrison format (2x VHH on C-terminus)
"""

from pathlib import Path
from typing import Optional
import yaml

from src.formatting.base import (
    BispecificConstruct,
    AntibodyChain,
    BispecificFormatter,
    SequenceLibrary,
)
from src.formatting.crossmab import CrossMabFormatter, assemble_crossmab
from src.formatting.fab_scfv import FabScFvFormatter, assemble_fab_scfv
from src.formatting.fab_vhh import FabVhhFormatter, assemble_fab_vhh
from src.formatting.igg_scfv import IggScfvFormatter, assemble_igg_scfv
from src.formatting.igg_vhh import IggVhhFormatter, assemble_igg_vhh
from src.utils.constants import parse_scfv, is_likely_scfv


# Registry of all formatters
FORMATTERS = {
    "crossmab": CrossMabFormatter,
    "fab_scfv": FabScFvFormatter,
    "fab_vhh": FabVhhFormatter,
    "igg_scfv": IggScfvFormatter,
    "igg_vhh": IggVhhFormatter,
}


def load_target_sequences(
    target_name: str,
    targets_file: Optional[str] = None,
) -> tuple[str, str, str]:
    """Load target arm sequences from placeholder_targets.yaml.

    Args:
        target_name: Name of target (e.g., 'trastuzumab', 'rituximab', 'cetuximab').
        targets_file: Path to targets YAML file. Defaults to data/frameworks/placeholder_targets.yaml.

    Returns:
        Tuple of (vh_sequence, vl_sequence, display_name).

    Raises:
        ValueError: If target_name not found in targets file.
        FileNotFoundError: If targets file not found.
    """
    if targets_file is None:
        # Try to find the file relative to common locations
        candidates = [
            Path("data/frameworks/placeholder_targets.yaml"),
            Path(__file__).parent.parent.parent / "data/frameworks/placeholder_targets.yaml",
        ]
        for candidate in candidates:
            if candidate.exists():
                targets_file = str(candidate)
                break

    if targets_file is None or not Path(targets_file).exists():
        raise FileNotFoundError(
            f"Targets file not found. Tried: {candidates if targets_file is None else targets_file}"
        )

    with open(targets_file, "r") as f:
        targets_data = yaml.safe_load(f)

    # Normalize target name for lookup
    target_key = target_name.lower().replace("-", "_").replace(" ", "_")

    if target_key not in targets_data:
        available = [k for k in targets_data.keys() if not k.startswith(("default", "usage"))]
        raise ValueError(
            f"Target '{target_name}' not found. Available targets: {available}"
        )

    target = targets_data[target_key]

    # Extract sequences (handle YAML multi-line strings)
    vh_seq = target.get("vh_sequence", "").replace("\n", "").replace(" ", "")
    vl_seq = target.get("vl_sequence", "").replace("\n", "").replace(" ", "")
    display_name = target.get("target", target.get("name", target_name.upper()))

    if not vh_seq or not vl_seq:
        raise ValueError(f"Target '{target_name}' missing vh_sequence or vl_sequence")

    return vh_seq, vl_seq, display_name


def get_formatter(format_type: str, sequence_library: SequenceLibrary = None) -> BispecificFormatter:
    """Get formatter instance by format type.

    Args:
        format_type: One of 'crossmab', 'fab_scfv', 'fab_vhh', 'igg_scfv', 'igg_vhh'.
        sequence_library: Custom sequence library.

    Returns:
        Formatter instance.

    Raises:
        ValueError: If format_type is not recognized.
    """
    if format_type not in FORMATTERS:
        raise ValueError(f"Unknown format: {format_type}. Available: {list(FORMATTERS.keys())}")

    return FORMATTERS[format_type](sequence_library)


def format_all(
    target_vh: str,
    target_vl: str,
    cd3_binder: str,
    cd3_binder_vl: str = None,
    name_prefix: str = "bispecific",
    target_name: str = "HER2",
    formats: list[str] = None,
    sequence_library: SequenceLibrary = None,
) -> dict[str, BispecificConstruct]:
    """Generate all bispecific formats for a CD3 binder.

    Args:
        target_vh: VH sequence for tumor target.
        target_vl: VL sequence for tumor target.
        cd3_binder: CD3 binder sequence (VHH, VH, or full scFv).
        cd3_binder_vl: VL for CD3 (required for scFv/CrossMab formats).
            If None but cd3_binder is a full scFv, will attempt to parse.
        name_prefix: Prefix for construct names.
        target_name: Name of tumor target.
        formats: List of formats to generate (default: all compatible).
        sequence_library: Custom sequence library.

    Returns:
        Dict mapping format names to BispecificConstruct instances.
    """
    results = {}
    library = sequence_library or SequenceLibrary()

    # If no VL provided, check if cd3_binder is a full scFv that we can parse
    effective_cd3_vh = cd3_binder
    effective_cd3_vl = cd3_binder_vl

    if effective_cd3_vl is None and is_likely_scfv(cd3_binder):
        parsed = parse_scfv(cd3_binder)
        if parsed:
            effective_cd3_vh, effective_cd3_vl = parsed
            print(f"  Parsed scFv into VH ({len(effective_cd3_vh)} aa) + VL ({len(effective_cd3_vl)} aa)")

    # Determine which formats are compatible
    is_vhh = effective_cd3_vl is None
    all_formats = formats or list(FORMATTERS.keys())

    for fmt in all_formats:
        # Skip formats that require VL if we only have VHH
        if is_vhh and fmt in ["crossmab", "fab_scfv", "igg_scfv"]:
            continue

        try:
            formatter = get_formatter(fmt, library)
            construct = formatter.assemble(
                target_vh=target_vh,
                target_vl=target_vl,
                cd3_binder=effective_cd3_vh,
                cd3_binder_vl=effective_cd3_vl,
                name=f"{name_prefix}_{fmt}",
                target_name=target_name,
            )
            results[fmt] = construct
        except Exception as e:
            # Log but don't fail on individual format errors
            print(f"Warning: Failed to generate {fmt} format: {e}")

    return results


__all__ = [
    "BispecificConstruct",
    "AntibodyChain",
    "BispecificFormatter",
    "SequenceLibrary",
    "CrossMabFormatter",
    "FabScFvFormatter",
    "FabVhhFormatter",
    "IggScfvFormatter",
    "IggVhhFormatter",
    "assemble_crossmab",
    "assemble_fab_scfv",
    "assemble_fab_vhh",
    "assemble_igg_scfv",
    "assemble_igg_vhh",
    "get_formatter",
    "format_all",
    "load_target_sequences",
    "FORMATTERS",
    "parse_scfv",
    "is_likely_scfv",
]
