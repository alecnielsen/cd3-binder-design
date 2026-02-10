"""Affinity proxy scoring using inverse folding models.

ProteinMPNN (MIT) and AntiFold (BSD-3) provide log-likelihood scores as
affinity proxies for de novo antibody designs. These are NOT direct affinity
predictors, but higher log-likelihoods correlate with better binding
(ProteinMPNN: Spearman r=0.27-0.41 on AbBiBench).

Usage:
    from src.analysis.affinity_scoring import batch_score_affinity
    results = batch_score_affinity("data/outputs/structures/cif/", design_ids, binder_types)
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import tempfile
import os


@dataclass
class AffinityResult:
    """Result from inverse folding affinity scoring."""

    design_id: str
    proteinmpnn_ll: Optional[float] = None  # NLL per residue (lower = better fit)
    antifold_ll: Optional[float] = None  # NLL per residue (lower = better fit)
    error: Optional[str] = None


# Cache loaded models to avoid reloading per candidate
_proteinmpnn_model = None
_antifold_model = None


def _cif_to_pdb(cif_path: str) -> str:
    """Convert CIF file to PDB format, returning path to temp PDB file.

    ProteinMPNN's parse_PDB only handles PDB format, not mmCIF.
    Uses BioPython's MMCIF2Dict + PDBIO for conversion.
    """
    from Bio.PDB import MMCIFParser, PDBIO
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("complex", cif_path)

    # Write to temp PDB file
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    pdb_io.save(tmp.name)
    return tmp.name


def _load_proteinmpnn_model():
    """Load ProteinMPNN model weights (cached)."""
    global _proteinmpnn_model
    if _proteinmpnn_model is not None:
        return _proteinmpnn_model

    import torch
    import pkg_resources
    from proteinmpnn.protein_mpnn_utils import ProteinMPNN

    device = torch.device("cpu")  # CPU is fine for scoring
    checkpoint_path = pkg_resources.resource_filename(
        "proteinmpnn", "data/vanilla_model_weights/v_48_020.pt"
    )
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)

    model = ProteinMPNN(
        ca_only=False,
        num_letters=21,
        node_features=128,
        edge_features=128,
        hidden_dim=128,
        num_encoder_layers=3,
        num_decoder_layers=3,
        k_neighbors=checkpoint["num_edges"],
    )
    model.to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    _proteinmpnn_model = (model, device)
    return _proteinmpnn_model


def score_proteinmpnn(cif_path: str, design_id: str) -> AffinityResult:
    """Score a complex using ProteinMPNN log-likelihood.

    ProteinMPNN scores how well the binder sequence "fits" the predicted
    complex structure. Lower negative log-likelihood = better structural fit.

    Args:
        cif_path: Path to CIF structure file from Boltz-2.
        design_id: Identifier for the design.

    Returns:
        AffinityResult with proteinmpnn_ll or error.
    """
    if not Path(cif_path).exists():
        return AffinityResult(design_id=design_id, error=f"CIF file not found: {cif_path}")

    try:
        from proteinmpnn.protein_mpnn_utils import parse_PDB, tied_featurize, _scores
    except ImportError:
        return AffinityResult(
            design_id=design_id,
            error="ProteinMPNN not installed. Install with: pip install proteinmpnn",
        )

    try:
        import torch

        model, device = _load_proteinmpnn_model()

        # Convert CIF to PDB if needed (ProteinMPNN only handles PDB format)
        parse_path = cif_path
        tmp_pdb = None
        if cif_path.endswith(".cif"):
            tmp_pdb = _cif_to_pdb(cif_path)
            parse_path = tmp_pdb

        try:
            pdb_dict_list = parse_PDB(parse_path)
        finally:
            if tmp_pdb and os.path.exists(tmp_pdb):
                os.unlink(tmp_pdb)

        if not pdb_dict_list:
            return AffinityResult(design_id=design_id, error="Failed to parse CIF file")

        pdb_dict = pdb_dict_list[0]

        # Identify chains — binder is chain B (or B+C for 3-chain), target is chain A
        all_chains = sorted(
            [k.replace("seq_chain_", "") for k in pdb_dict.keys() if k.startswith("seq_chain_")]
        )

        if len(all_chains) < 2:
            return AffinityResult(design_id=design_id, error=f"Need >=2 chains, found {len(all_chains)}")

        # Score binder chain(s) conditioned on the full complex
        # Chain A = target (fixed), B (and C if present) = binder (score these)
        target_chains = [all_chains[0]]
        binder_chains = all_chains[1:]

        # Build chain_id_dict: {name: (chains_to_design, chains_to_keep_fixed)}
        chain_id_dict = {pdb_dict["name"]: (binder_chains, target_chains)}

        # Featurize
        X, S, mask, lengths, chain_M, chain_encoding_all, chain_list_list, \
            visible_list_list, masked_list_list, masked_chain_length_list_list, \
            chain_M_pos, omit_AA_mask, residue_idx, dihedral_mask, tied_pos_list_of_lists, \
            pssm_coef, pssm_bias, pssm_log_odds_all, bias_by_res, tied_beta = \
            tied_featurize(
                [pdb_dict], device, chain_id_dict,
                fixed_position_dict=None, omit_AA_dict=None,
                tied_positions_dict=None, pssm_dict=None,
                bias_by_res_dict=None, ca_only=False,
            )

        # Forward pass
        with torch.no_grad():
            randn = torch.randn(chain_M.shape, device=device)
            log_probs = model(
                X, S, mask, chain_M, residue_idx, chain_encoding_all, randn,
            )
            # Score only the binder chain(s)
            mask_for_loss = mask * chain_M
            scores = _scores(S, log_probs, mask_for_loss)

        ll = float(scores[0].cpu().numpy())
        return AffinityResult(design_id=design_id, proteinmpnn_ll=ll)

    except Exception as e:
        return AffinityResult(design_id=design_id, error=f"ProteinMPNN error: {e}")


def score_antifold(cif_path: str, design_id: str, binder_type: str = "vhh") -> AffinityResult:
    """Score a complex using AntiFold log-likelihood.

    AntiFold is an antibody-specific inverse folding model trained on
    paired antibody structures. It supports nanobodies via chain specification.

    Args:
        cif_path: Path to CIF structure file.
        design_id: Identifier for the design.
        binder_type: "vhh" for nanobody, "scfv"/"fab" for paired.

    Returns:
        AffinityResult with antifold_ll or error.
    """
    if not Path(cif_path).exists():
        return AffinityResult(design_id=design_id, error=f"CIF file not found: {cif_path}")

    try:
        from antifold.antiscripts import load_model, get_pdbs_logits
    except ImportError:
        return AffinityResult(
            design_id=design_id,
            error="AntiFold not installed. Install with: pip install antifold",
        )

    try:
        import pandas as pd

        global _antifold_model
        if _antifold_model is None:
            _antifold_model = load_model(checkpoint_path="")

        cif_file = Path(cif_path)
        pdb_name = cif_file.stem  # e.g., "fab_1XIW_0008"

        # AntiFold needs a CSV mapping PDB names to chain IDs
        # Boltz-2 CIF: chain A = target, chain B = binder (or B=VH, C=VL)
        # Boltz-2 CIF chain layout:
        #   2-chain (VHH/scFv): A=target, B=binder
        #   3-chain (Fab): A=target, B=VH, C=VL
        is_3chain = "_3chain" in str(cif_path)

        if binder_type == "vhh" or (binder_type in ("scfv", "fab") and not is_3chain):
            # Single binder chain B — score as nanobody/single-chain
            # custom_chain_mode=True required: skips Hchain+Lchain column validation
            df_csv = pd.DataFrame({
                "pdb": [pdb_name],
                "Hchain": ["B"],
                "Lchain": [None],
            })
            use_custom_chain = True
        else:
            # 3-chain: B=VH, C=VL
            df_csv = pd.DataFrame({
                "pdb": [pdb_name],
                "Hchain": ["B"],
                "Lchain": ["C"],
            })
            use_custom_chain = False

        # Run inference
        df_logits_list = get_pdbs_logits(
            _antifold_model,
            pdbs_csv_or_dataframe=df_csv,
            pdb_dir=str(cif_file.parent),
            batch_size=1,
            custom_chain_mode=use_custom_chain,
            num_threads=0,
            save_flag=False,
            seed=42,
        )

        if not df_logits_list or len(df_logits_list) == 0:
            return AffinityResult(design_id=design_id, error="AntiFold returned no results")

        df_logits = df_logits_list[0]

        # Compute global score from logits
        # The logits DataFrame has amino acid columns — compute NLL from these
        import torch
        import torch.nn.functional as F

        aa_cols = list("ACDEFGHIKLMNPQRSTVWY")
        available_cols = [c for c in aa_cols if c in df_logits.columns]
        if not available_cols:
            return AffinityResult(design_id=design_id, error="No amino acid columns in AntiFold output")

        logits_tensor = torch.tensor(df_logits[available_cols].values, dtype=torch.float32)
        log_probs = F.log_softmax(logits_tensor, dim=-1)

        # Get the actual residues — column name is "pdb_res" in AntiFold output
        res_col = "pdb_res" if "pdb_res" in df_logits.columns else "wt"
        if res_col not in df_logits.columns:
            return AffinityResult(design_id=design_id, error="No residue column in AntiFold logits")
        wt_residues = df_logits[res_col].tolist()

        # Compute NLL for the actual sequence
        total_nll = 0.0
        count = 0
        aa_to_idx = {aa: i for i, aa in enumerate(available_cols)}
        for i, res in enumerate(wt_residues):
            if res in aa_to_idx:
                idx = aa_to_idx[res]
                total_nll -= float(log_probs[i, idx])
                count += 1

        if count == 0:
            return AffinityResult(design_id=design_id, error="No valid residues for scoring")

        mean_nll = total_nll / count
        return AffinityResult(design_id=design_id, antifold_ll=mean_nll)

    except Exception as e:
        return AffinityResult(design_id=design_id, error=f"AntiFold error: {e}")


def batch_score_affinity(
    cif_dir: str,
    design_ids: list,
    binder_types: list,
    run_proteinmpnn: bool = True,
    run_antifold: bool = True,
) -> list:
    """Batch affinity scoring for multiple designs.

    Runs both ProteinMPNN and AntiFold on each design, merging results.
    Graceful error handling: returns error message if a tool is not installed.

    Args:
        cif_dir: Directory containing CIF files (named {design_id}.cif).
        design_ids: List of design identifiers.
        binder_types: Corresponding binder types ("vhh" or "scfv").
        run_proteinmpnn: Whether to run ProteinMPNN scoring.
        run_antifold: Whether to run AntiFold scoring.

    Returns:
        List of AffinityResult objects with merged scores.
    """
    cif_base = Path(cif_dir)
    results = []

    for design_id, binder_type in zip(design_ids, binder_types):
        cif_path = cif_base / f"{design_id}.cif"
        merged = AffinityResult(design_id=design_id)
        errors = []

        if run_proteinmpnn:
            mpnn_result = score_proteinmpnn(str(cif_path), design_id)
            if mpnn_result.proteinmpnn_ll is not None:
                merged.proteinmpnn_ll = mpnn_result.proteinmpnn_ll
            elif mpnn_result.error:
                errors.append(mpnn_result.error)

        if run_antifold:
            af_result = score_antifold(str(cif_path), design_id, binder_type)
            if af_result.antifold_ll is not None:
                merged.antifold_ll = af_result.antifold_ll
            elif af_result.error:
                errors.append(af_result.error)

        if errors and merged.proteinmpnn_ll is None and merged.antifold_ll is None:
            merged.error = "; ".join(errors)

        results.append(merged)

    return results
