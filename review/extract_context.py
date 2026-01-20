#!/usr/bin/env python3
"""
Extract context slices from README.md for scientific code review.

This script reads context_slices.yaml and extracts the specified line ranges
from README.md, combining them into a single context string for a given module.

Usage:
    python extract_context.py <module_name>
    python extract_context.py analysis
    python extract_context.py --list  # List available modules
"""

import argparse
import sys
from pathlib import Path

import yaml


def load_context_slices(yaml_path: Path) -> dict:
    """Load context slice definitions from YAML."""
    with open(yaml_path, "r") as f:
        return yaml.safe_load(f)


def extract_lines(readme_path: Path, start: int, end: int) -> list[str]:
    """Extract lines from README (1-indexed, inclusive)."""
    with open(readme_path, "r") as f:
        lines = f.readlines()
    # Convert to 0-indexed
    return lines[start - 1 : end]


def get_module_context(
    module_name: str, slices: dict, readme_path: Path
) -> tuple[str, list[str]]:
    """
    Extract combined context for a module.

    Returns:
        tuple: (combined_context_string, list_of_source_paths)
    """
    if module_name not in slices:
        available = [k for k in slices.keys() if not k.startswith("_")]
        raise ValueError(
            f"Unknown module: {module_name}. Available: {', '.join(available)}"
        )

    module_config = slices[module_name]
    sections = module_config.get("readme_sections", [])
    source_paths = module_config.get("source_paths", [])
    description = module_config.get("description", "")

    context_parts = []

    # Add module description header
    context_parts.append(f"# README Context for Module: {module_name}")
    context_parts.append(f"# Description: {description}")
    context_parts.append("")

    # Extract each section
    for section in sections:
        name = section["name"]
        start = section["start"]
        end = section["end"]

        context_parts.append(f"## {name} (lines {start}-{end})")
        context_parts.append("")

        lines = extract_lines(readme_path, start, end)
        context_parts.extend(line.rstrip() for line in lines)
        context_parts.append("")
        context_parts.append("---")
        context_parts.append("")

    return "\n".join(context_parts), source_paths


def list_modules(slices: dict) -> None:
    """Print available modules and their descriptions."""
    print("Available modules for review:\n")
    for name, config in slices.items():
        if isinstance(config, dict) and "description" in config:
            desc = config.get("description", "No description")
            sections = config.get("readme_sections", [])
            sources = config.get("source_paths", [])
            print(f"  {name}:")
            print(f"    Description: {desc}")
            print(f"    README sections: {len(sections)}")
            print(f"    Source paths: {', '.join(sources)}")
            print()


def main():
    parser = argparse.ArgumentParser(
        description="Extract README context slices for module review"
    )
    parser.add_argument(
        "module",
        nargs="?",
        help="Module name to extract context for (e.g., analysis, formatting)",
    )
    parser.add_argument(
        "--list", "-l", action="store_true", help="List available modules"
    )
    parser.add_argument(
        "--slices",
        default=None,
        help="Path to context_slices.yaml (default: same directory as script)",
    )
    parser.add_argument(
        "--readme",
        default=None,
        help="Path to README.md (default: parent directory of script)",
    )
    parser.add_argument(
        "--sources-only",
        action="store_true",
        help="Only output source paths, not README context",
    )

    args = parser.parse_args()

    # Resolve paths relative to script location
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent

    slices_path = Path(args.slices) if args.slices else script_dir / "context_slices.yaml"
    readme_path = Path(args.readme) if args.readme else repo_root / "README.md"

    if not slices_path.exists():
        print(f"Error: Context slices file not found: {slices_path}", file=sys.stderr)
        sys.exit(1)

    if not readme_path.exists():
        print(f"Error: README.md not found: {readme_path}", file=sys.stderr)
        sys.exit(1)

    slices = load_context_slices(slices_path)

    if args.list:
        list_modules(slices)
        sys.exit(0)

    if not args.module:
        parser.print_help()
        sys.exit(1)

    try:
        context, source_paths = get_module_context(args.module, slices, readme_path)

        if args.sources_only:
            for path in source_paths:
                print(path)
        else:
            print(context)

    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
