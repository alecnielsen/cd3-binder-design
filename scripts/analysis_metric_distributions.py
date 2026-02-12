#!/usr/bin/env python3
"""Statistical analysis of 192-design metric distributions and sequence diversity.

Analyzes the 100x scale run to understand:
1. Distribution of composite scores, interface areas, contacts, humanness
2. How the top 10 differ from the overall population
3. Sequence diversity among final candidates (identity matrix, CDR-H3)
4. Scaling diminishing returns (10x vs 100x comparison)
"""

import json
import sys
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(__file__).resolve().parent.parent
STRUCTURES = BASE / "data/outputs/structures/candidates_with_structures.json"
FILTERED = BASE / "data/outputs/filtered/filtered_candidates.json"
DENOVO_100X = BASE / "data/outputs/denovo/denovo_results_20260205_173718.json"
DENOVO_10X = BASE / "data/outputs/denovo/denovo_results_20260205_101836.json"
OUT_DIR = BASE / "data/outputs/analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_all_candidates():
    """Load the full 192-design dataset with structure predictions."""
    with open(STRUCTURES) as f:
        data = json.load(f)
    return data["candidates"]


def load_filtered():
    """Load the 10 final candidates."""
    with open(FILTERED) as f:
        data = json.load(f)
    return data["candidates"]


def load_denovo(path):
    """Load raw denovo design results."""
    with open(path) as f:
        data = json.load(f)
    designs = []
    for key in ("vhh_designs", "fab_designs"):
        if key in data:
            designs.extend(data[key])
    return designs


def build_dataframe(candidates):
    """Convert candidates list to a tidy DataFrame."""
    rows = []
    for c in candidates:
        sp = c.get("structure_prediction") or {}
        row = {
            "design_id": c.get("design_id", ""),
            "binder_type": c.get("binder_type", ""),
            "sequence": c.get("sequence", ""),
            "seq_len": len(c.get("sequence", "")),
            "iptm": c.get("iptm", np.nan),
            "ptm": sp.get("ptm", np.nan),
            "plddt": sp.get("plddt_mean", np.nan),
            "interface_area": sp.get("interface_area", np.nan),
            "num_contacts": sp.get("num_contacts", np.nan),
            "has_structure": bool(sp),
        }
        rows.append(row)
    df = pd.DataFrame(rows)
    return df


def build_filtered_df(candidates):
    """Build DataFrame from filtered candidates with all scored fields."""
    rows = []
    for c in candidates:
        row = {
            "design_id": c.get("candidate_id", ""),
            "binder_type": c.get("binder_type", ""),
            "sequence": c.get("sequence", ""),
            "composite_score": c.get("composite_score", np.nan),
            "rank": c.get("rank", np.nan),
            "interface_area": c["binding"]["interface_area"],
            "num_contacts": c["binding"]["num_contacts"],
            "humanness": c["humanness"]["oasis_score_mean"],
            "net_charge": c["developability"]["net_charge"],
            "pi": c["developability"]["isoelectric_point"],
            "epitope_class": c["epitope"]["epitope_class"],
            "okt3_overlap": c["epitope"]["okt3_overlap"],
            "n_risk_flags": len(c.get("risk_flags", [])),
        }
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Sequence analysis helpers
# ---------------------------------------------------------------------------
def pairwise_identity(seq1, seq2):
    """Compute pairwise sequence identity (shorter / longer alignment)."""
    # Simple global: count matches at each position up to min length
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    if max_len == 0:
        return 0.0
    matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
    return matches / max_len


def extract_cdrh3_approx(sequence):
    """Extract approximate CDR-H3 from a VH/VHH sequence.

    Uses conserved motifs flanking CDR-H3. Tries multiple patterns
    since VHH and Fab frameworks differ.
    """
    import re
    # Standard Fab: C[A/V][R/K]...WG[Q/K]G
    match = re.search(r"C[AVST][RK](.{3,30}?)W[GD][QK]G", sequence)
    if match:
        return match.group(1)
    # VHH/broader: YY[C/W/E]...WG  (VHH often has YYE, YYA patterns)
    match = re.search(r"YY[CWEA](.{3,30}?)W[GD][QK]G", sequence)
    if match:
        return match.group(1)
    # Ultra-broad fallback: YY...WG
    match = re.search(r"YY(.{3,30}?)WG", sequence)
    if match:
        return match.group(1)
    return None


# ---------------------------------------------------------------------------
# Part 1: Metric distributions (all 192 designs)
# ---------------------------------------------------------------------------
def analyze_distributions(df_all, df_top10):
    """Statistical summary and distribution plots for all designs."""
    # Only designs with successful structure prediction
    df_pred = df_all[df_all["has_structure"]].copy()

    print("=" * 70)
    print("PART 1: METRIC DISTRIBUTIONS (192 designs, {} with structures)".format(len(df_pred)))
    print("=" * 70)

    metrics = {
        "iptm": "ipTM (design-time)",
        "interface_area": "Interface Area (A^2)",
        "num_contacts": "Residue Contacts",
        "ptm": "pTM (Boltz-2)",
        "plddt": "pLDDT mean",
    }

    # Summary stats
    print("\n--- Summary Statistics (all predicted designs) ---")
    for col, label in metrics.items():
        vals = df_pred[col].dropna()
        if len(vals) == 0:
            continue
        print(f"\n{label}:")
        print(f"  n={len(vals)}  mean={vals.mean():.3f}  std={vals.std():.3f}")
        print(f"  min={vals.min():.3f}  25%={vals.quantile(0.25):.3f}  "
              f"median={vals.median():.3f}  75%={vals.quantile(0.75):.3f}  max={vals.max():.3f}")

    # By type
    print("\n--- By Binder Type ---")
    for btype in ["vhh", "fab"]:
        sub = df_pred[df_pred["binder_type"] == btype]
        print(f"\n{btype.upper()} (n={len(sub)}):")
        for col in ["iptm", "interface_area", "num_contacts"]:
            vals = sub[col].dropna()
            if len(vals) > 0:
                print(f"  {col}: mean={vals.mean():.3f}  std={vals.std():.3f}  "
                      f"median={vals.median():.3f}  max={vals.max():.3f}")

    # Top 10 vs population comparison
    print("\n--- Top 10 vs Population ---")
    top10_ids = set(df_top10["design_id"])
    df_pred["is_top10"] = df_pred["design_id"].isin(top10_ids)

    for col in ["interface_area", "num_contacts", "iptm"]:
        pop = df_pred[~df_pred["is_top10"]][col].dropna()
        top = df_pred[df_pred["is_top10"]][col].dropna()
        if len(pop) > 0 and len(top) > 0:
            # Mann-Whitney U test (non-parametric)
            stat, pval = stats.mannwhitneyu(top, pop, alternative="greater")
            effect = top.mean() - pop.mean()
            print(f"\n{col}:")
            print(f"  Population: mean={pop.mean():.1f}  Top10: mean={top.mean():.1f}")
            print(f"  Difference: +{effect:.1f}  (Mann-Whitney p={pval:.4f})")
            # Percentile of top 10 within population
            percentiles = [stats.percentileofscore(pop, v) for v in top]
            print(f"  Top10 percentiles: {[f'{p:.0f}' for p in sorted(percentiles)]}")

    return df_pred


def plot_distributions(df_pred, df_top10):
    """Generate distribution plots."""
    top10_ids = set(df_top10["design_id"])
    df_pred["is_top10"] = df_pred["design_id"].isin(top10_ids)

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle("CD3 Binder Design - Metric Distributions (100x Scale Run)", fontsize=14, y=0.98)

    # 1. ipTM distribution
    ax = axes[0, 0]
    for btype, color in [("vhh", "#2196F3"), ("fab", "#FF9800")]:
        vals = df_pred[df_pred["binder_type"] == btype]["iptm"].dropna()
        ax.hist(vals, bins=20, alpha=0.6, color=color, label=btype.upper(), edgecolor="white")
    # Mark top 10
    top10_iptm = df_pred[df_pred["is_top10"]]["iptm"].dropna()
    for v in top10_iptm:
        ax.axvline(v, color="red", alpha=0.5, linewidth=1, linestyle="--")
    ax.set_xlabel("ipTM (design-time)")
    ax.set_ylabel("Count")
    ax.set_title("ipTM Distribution")
    ax.legend()

    # 2. Interface area distribution
    ax = axes[0, 1]
    for btype, color in [("vhh", "#2196F3"), ("fab", "#FF9800")]:
        vals = df_pred[df_pred["binder_type"] == btype]["interface_area"].dropna()
        ax.hist(vals, bins=20, alpha=0.6, color=color, label=btype.upper(), edgecolor="white")
    top10_ia = df_pred[df_pred["is_top10"]]["interface_area"].dropna()
    for v in top10_ia:
        ax.axvline(v, color="red", alpha=0.5, linewidth=1, linestyle="--")
    ax.axvline(2060, color="black", linewidth=2, linestyle=":", label="Threshold (2060)")
    ax.set_xlabel("Interface Area (A^2)")
    ax.set_ylabel("Count")
    ax.set_title("Interface Area Distribution")
    ax.legend(fontsize=8)

    # 3. Contact count distribution
    ax = axes[0, 2]
    for btype, color in [("vhh", "#2196F3"), ("fab", "#FF9800")]:
        vals = df_pred[df_pred["binder_type"] == btype]["num_contacts"].dropna()
        ax.hist(vals, bins=20, alpha=0.6, color=color, label=btype.upper(), edgecolor="white")
    top10_nc = df_pred[df_pred["is_top10"]]["num_contacts"].dropna()
    for v in top10_nc:
        ax.axvline(v, color="red", alpha=0.5, linewidth=1, linestyle="--")
    ax.axvline(28, color="black", linewidth=2, linestyle=":", label="Threshold (28)")
    ax.set_xlabel("Residue Contacts")
    ax.set_ylabel("Count")
    ax.set_title("Contact Count Distribution")
    ax.legend(fontsize=8)

    # 4. Interface area vs contacts scatter
    ax = axes[1, 0]
    df_pop = df_pred[~df_pred["is_top10"]]
    df_top = df_pred[df_pred["is_top10"]]
    for btype, marker in [("vhh", "o"), ("fab", "s")]:
        sub = df_pop[df_pop["binder_type"] == btype]
        ax.scatter(sub["interface_area"], sub["num_contacts"], alpha=0.3,
                   marker=marker, s=20, label=f"{btype.upper()} (pop)")
    for btype, marker in [("vhh", "o"), ("fab", "s")]:
        sub = df_top[df_top["binder_type"] == btype]
        ax.scatter(sub["interface_area"], sub["num_contacts"], alpha=0.9,
                   marker=marker, s=80, edgecolors="red", linewidths=2,
                   label=f"{btype.upper()} (top10)")
    ax.axvline(2060, color="black", linestyle=":", alpha=0.5)
    ax.axhline(28, color="black", linestyle=":", alpha=0.5)
    ax.set_xlabel("Interface Area (A^2)")
    ax.set_ylabel("Residue Contacts")
    ax.set_title("Interface Area vs Contacts")
    ax.legend(fontsize=7)

    # 5. pTM vs pLDDT scatter
    ax = axes[1, 1]
    for btype, marker, color in [("vhh", "o", "#2196F3"), ("fab", "s", "#FF9800")]:
        sub = df_pred[df_pred["binder_type"] == btype]
        ax.scatter(sub["ptm"], sub["plddt"], alpha=0.4, marker=marker,
                   s=20, color=color, label=btype.upper())
    ax.set_xlabel("pTM")
    ax.set_ylabel("pLDDT mean")
    ax.set_title("Structural Confidence")
    ax.legend()

    # 6. Composite score (top 10 bar chart)
    ax = axes[1, 2]
    df_top10_sorted = df_top10.sort_values("composite_score", ascending=True)
    colors = ["#FF9800" if t == "fab" else "#2196F3"
              for t in df_top10_sorted["binder_type"]]
    bars = ax.barh(range(len(df_top10_sorted)), df_top10_sorted["composite_score"],
                   color=colors, edgecolor="white")
    ax.set_yticks(range(len(df_top10_sorted)))
    ax.set_yticklabels(df_top10_sorted["design_id"], fontsize=8)
    ax.set_xlabel("Composite Score")
    ax.set_title("Top 10 Candidates")
    # Legend
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#FF9800", label="Fab"),
                       Patch(color="#2196F3", label="VHH")], fontsize=8)

    plt.tight_layout()
    path = OUT_DIR / "metric_distributions.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {path}")
    return path


# ---------------------------------------------------------------------------
# Part 2: Sequence diversity
# ---------------------------------------------------------------------------
def analyze_diversity(df_top10):
    """Sequence diversity analysis of top 10 candidates."""
    print("\n" + "=" * 70)
    print("PART 2: SEQUENCE DIVERSITY (Top 10 Candidates)")
    print("=" * 70)

    seqs = df_top10["sequence"].tolist()
    ids = df_top10["design_id"].tolist()
    n = len(seqs)

    # Pairwise identity matrix
    identity_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                identity_matrix[i, j] = 1.0
            else:
                identity_matrix[i, j] = pairwise_identity(seqs[i], seqs[j])

    print("\n--- Pairwise Sequence Identity Matrix ---")
    id_df = pd.DataFrame(identity_matrix, index=ids, columns=ids)
    # Print condensed
    for i in range(n):
        vals = [f"{identity_matrix[i, j]:.2f}" for j in range(n)]
        short_id = ids[i].replace("fab_", "F:").replace("vhh_", "V:")
        print(f"  {short_id:18s} " + " ".join(vals))

    # Summary stats on pairwise identities (off-diagonal)
    off_diag = []
    for i, j in combinations(range(n), 2):
        off_diag.append(identity_matrix[i, j])
    off_diag = np.array(off_diag)

    print(f"\n--- Pairwise Identity Summary ---")
    print(f"  Mean: {off_diag.mean():.3f}")
    print(f"  Std:  {off_diag.std():.3f}")
    print(f"  Min:  {off_diag.min():.3f}  (most diverse pair)")
    print(f"  Max:  {off_diag.max():.3f}  (most similar pair)")

    # Same-type vs cross-type
    types = df_top10["binder_type"].tolist()
    same_type, cross_type = [], []
    for i, j in combinations(range(n), 2):
        if types[i] == types[j]:
            same_type.append(identity_matrix[i, j])
        else:
            cross_type.append(identity_matrix[i, j])

    if same_type:
        print(f"\n  Same-type pairs:  mean={np.mean(same_type):.3f}  (n={len(same_type)})")
    if cross_type:
        print(f"  Cross-type pairs: mean={np.mean(cross_type):.3f}  (n={len(cross_type)})")

    # CDR-H3 analysis
    print("\n--- CDR-H3 Analysis ---")
    cdr_h3s = []
    for i, seq in enumerate(seqs):
        cdr = extract_cdrh3_approx(seq)
        cdr_h3s.append(cdr)
        if cdr:
            print(f"  {ids[i]:22s}  CDR-H3: {cdr}  (len={len(cdr)})")
        else:
            print(f"  {ids[i]:22s}  CDR-H3: not found")

    # CDR-H3 pairwise identity
    valid_cdrs = [(i, cdr_h3s[i]) for i in range(n) if cdr_h3s[i] is not None]
    if len(valid_cdrs) >= 2:
        print(f"\n--- CDR-H3 Pairwise Identity ---")
        cdr_ids = []
        for (i, ci), (j, cj) in combinations(valid_cdrs, 2):
            ident = pairwise_identity(ci, cj)
            cdr_ids.append(ident)
            if ident > 0.6:
                print(f"  {ids[i]} vs {ids[j]}: {ident:.2f} ***")
        print(f"\n  CDR-H3 mean identity: {np.mean(cdr_ids):.3f}")
        print(f"  CDR-H3 max identity:  {np.max(cdr_ids):.3f}")
        print(f"  CDR-H3 unique count:  {len(set(c for _, c in valid_cdrs))}/{len(valid_cdrs)}")

    # Sequence length stats
    print(f"\n--- Sequence Lengths ---")
    for btype in ["fab", "vhh"]:
        sub = df_top10[df_top10["binder_type"] == btype]
        if len(sub) > 0:
            lens = sub["sequence"].str.len()
            print(f"  {btype.upper()}: mean={lens.mean():.0f}  range=[{lens.min()}, {lens.max()}]")

    return identity_matrix, id_df


def plot_diversity(df_top10, identity_matrix):
    """Generate diversity plots."""
    ids = df_top10["design_id"].tolist()
    n = len(ids)
    short_ids = [i.replace("fab_", "F:").replace("vhh_", "V:") for i in ids]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("CD3 Binder Design - Sequence Diversity Analysis", fontsize=14, y=1.02)

    # 1. Identity heatmap
    ax = axes[0]
    im = ax.imshow(identity_matrix, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(n))
    ax.set_xticklabels(short_ids, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n))
    ax.set_yticklabels(short_ids, fontsize=7)
    ax.set_title("Pairwise Sequence Identity")
    # Add text annotations
    for i in range(n):
        for j in range(n):
            color = "white" if identity_matrix[i, j] > 0.6 else "black"
            ax.text(j, i, f"{identity_matrix[i, j]:.2f}", ha="center", va="center",
                    fontsize=6, color=color)
    plt.colorbar(im, ax=ax, shrink=0.8)

    # 2. Hierarchical clustering dendrogram
    ax = axes[1]
    # Convert identity to distance
    dist_matrix = 1.0 - identity_matrix
    np.fill_diagonal(dist_matrix, 0)
    # Make symmetric (handle float precision)
    dist_matrix = (dist_matrix + dist_matrix.T) / 2
    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method="average")
    dendrogram(Z, labels=short_ids, ax=ax, leaf_rotation=45, leaf_font_size=8)
    ax.set_title("Hierarchical Clustering (avg linkage)")
    ax.set_ylabel("Sequence Distance (1 - identity)")

    # Determine clusters at 0.7 distance (30% identity)
    clusters = fcluster(Z, t=0.7, criterion="distance")
    print(f"\n--- Clustering (distance threshold 0.7) ---")
    for cl in sorted(set(clusters)):
        members = [ids[i] for i in range(n) if clusters[i] == cl]
        print(f"  Cluster {cl}: {members}")

    # 3. CDR-H3 length vs humanness
    ax = axes[2]
    cdr_lengths = []
    for seq in df_top10["sequence"]:
        cdr = extract_cdrh3_approx(seq)
        cdr_lengths.append(len(cdr) if cdr else 0)
    df_top10 = df_top10.copy()
    df_top10["cdr_h3_len"] = cdr_lengths

    for btype, color, marker in [("fab", "#FF9800", "s"), ("vhh", "#2196F3", "o")]:
        sub = df_top10[df_top10["binder_type"] == btype]
        ax.scatter(sub["cdr_h3_len"], sub["humanness"], s=100, color=color,
                   marker=marker, edgecolors="black", linewidths=1, label=btype.upper(), zorder=3)
        for _, row in sub.iterrows():
            short = row["design_id"].replace("fab_", "F:").replace("vhh_", "V:")
            ax.annotate(short, (row["cdr_h3_len"], row["humanness"]),
                        textcoords="offset points", xytext=(5, 5), fontsize=6)
    ax.set_xlabel("CDR-H3 Length (aa)")
    ax.set_ylabel("Humanness (Sapiens)")
    ax.set_title("CDR-H3 Length vs Humanness")
    ax.axhline(0.8, color="gray", linestyle=":", alpha=0.5, label="Threshold (0.8)")
    ax.legend(fontsize=8)

    plt.tight_layout()
    path = OUT_DIR / "sequence_diversity.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {path}")
    return path


# ---------------------------------------------------------------------------
# Part 3: Scaling analysis (10x vs 100x)
# ---------------------------------------------------------------------------
def analyze_scaling(df_all):
    """Compare what we'd get from 10x vs 100x sampling."""
    print("\n" + "=" * 70)
    print("PART 3: SCALING / DIMINISHING RETURNS ANALYSIS")
    print("=" * 70)

    df_pred = df_all[df_all["has_structure"]].copy()

    # Simulate picking top-N by interface_area (primary filter)
    # Sort by interface_area descending
    df_sorted = df_pred.sort_values("interface_area", ascending=False)

    # What would sampling N designs yield in terms of best interface_area?
    sample_sizes = [10, 20, 50, 100, 150, len(df_sorted)]
    print("\n--- Best Interface Area by Sample Size ---")
    for n in sample_sizes:
        if n > len(df_sorted):
            continue
        top_n = df_sorted.head(n)
        print(f"  Top {n:3d}: best={top_n['interface_area'].max():.0f}  "
              f"mean_top10={top_n.head(10)['interface_area'].mean():.0f}  "
              f"passing_threshold={len(top_n[top_n['interface_area'] >= 2060])}")

    # Bootstrap analysis: resample to estimate expected improvement at larger N
    print("\n--- Bootstrap: Expected Top-1 Interface Area at Different Scales ---")
    rng = np.random.default_rng(42)
    all_vals = df_pred["interface_area"].dropna().values
    scales = [50, 100, 200, 500, 1000]
    n_boot = 1000

    print(f"  {'Scale':>6s}  {'Mean Top-1':>10s}  {'95% CI':>18s}  {'Mean Top-10':>11s}")
    for scale in scales:
        top1_vals = []
        top10_means = []
        for _ in range(n_boot):
            sample = rng.choice(all_vals, size=min(scale, len(all_vals)), replace=True)
            sample_sorted = np.sort(sample)[::-1]
            top1_vals.append(sample_sorted[0])
            top10_means.append(sample_sorted[:min(10, len(sample_sorted))].mean())
        ci_lo, ci_hi = np.percentile(top1_vals, [2.5, 97.5])
        print(f"  {scale:6d}  {np.mean(top1_vals):10.0f}  [{ci_lo:7.0f}, {ci_hi:7.0f}]  "
              f"{np.mean(top10_means):11.0f}")

    # Same for contacts
    print("\n--- Bootstrap: Expected Top-1 Contacts at Different Scales ---")
    all_contacts = df_pred["num_contacts"].dropna().values
    print(f"  {'Scale':>6s}  {'Mean Top-1':>10s}  {'95% CI':>18s}")
    for scale in scales:
        top1_vals = []
        for _ in range(n_boot):
            sample = rng.choice(all_contacts, size=min(scale, len(all_contacts)), replace=True)
            top1_vals.append(np.max(sample))
        ci_lo, ci_hi = np.percentile(top1_vals, [2.5, 97.5])
        print(f"  {scale:6d}  {np.mean(top1_vals):10.0f}  [{ci_lo:7.0f}, {ci_hi:7.0f}]")

    return df_sorted


def plot_scaling(df_all):
    """Scaling analysis plots."""
    df_pred = df_all[df_all["has_structure"]].copy()

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("CD3 Binder Design - Scaling Analysis", fontsize=14, y=1.02)

    # 1. Cumulative max interface area
    ax = axes[0]
    for btype, color in [("vhh", "#2196F3"), ("fab", "#FF9800"), ("all", "#4CAF50")]:
        if btype == "all":
            vals = df_pred["interface_area"].dropna().values
        else:
            vals = df_pred[df_pred["binder_type"] == btype]["interface_area"].dropna().values
        rng = np.random.default_rng(42)
        # Show cumulative max as we sample more designs
        n_boot = 100
        x = np.arange(1, len(vals) + 1)
        cum_maxes = np.zeros((n_boot, len(vals)))
        for b in range(n_boot):
            perm = rng.permutation(vals)
            cum_maxes[b] = np.maximum.accumulate(perm)
        mean_line = cum_maxes.mean(axis=0)
        lo = np.percentile(cum_maxes, 5, axis=0)
        hi = np.percentile(cum_maxes, 95, axis=0)
        label = btype.upper() if btype != "all" else "Combined"
        ax.plot(x, mean_line, color=color, label=label, linewidth=2)
        ax.fill_between(x, lo, hi, color=color, alpha=0.15)
    ax.axhline(2060, color="black", linestyle=":", alpha=0.5, label="Filter threshold")
    ax.set_xlabel("Number of Designs Sampled")
    ax.set_ylabel("Best Interface Area (A^2)")
    ax.set_title("Cumulative Max Interface Area")
    ax.legend(fontsize=8)

    # 2. Pass rate at different thresholds
    ax = axes[1]
    thresholds_ia = np.arange(500, 5500, 200)
    for btype, color in [("vhh", "#2196F3"), ("fab", "#FF9800")]:
        vals = df_pred[df_pred["binder_type"] == btype]["interface_area"].dropna().values
        rates = [(vals >= t).sum() / len(vals) * 100 for t in thresholds_ia]
        ax.plot(thresholds_ia, rates, color=color, label=btype.upper(), linewidth=2)
    ax.axvline(2060, color="black", linestyle=":", alpha=0.5)
    ax.set_xlabel("Interface Area Threshold (A^2)")
    ax.set_ylabel("Pass Rate (%)")
    ax.set_title("Filter Pass Rate vs Threshold")
    ax.legend()

    # 3. Score distribution by type (violin)
    ax = axes[2]
    data_vhh = df_pred[df_pred["binder_type"] == "vhh"]["interface_area"].dropna()
    data_fab = df_pred[df_pred["binder_type"] == "fab"]["interface_area"].dropna()
    parts = ax.violinplot([data_vhh, data_fab], positions=[1, 2], showmeans=True, showmedians=True)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(["#2196F3", "#FF9800"][i])
        pc.set_alpha(0.6)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["VHH", "Fab"])
    ax.set_ylabel("Interface Area (A^2)")
    ax.set_title("VHH vs Fab Interface Distribution")
    ax.axhline(2060, color="black", linestyle=":", alpha=0.5, label="Threshold")
    ax.legend()

    plt.tight_layout()
    path = OUT_DIR / "scaling_analysis.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {path}")
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Loading data...")
    all_candidates = load_all_candidates()
    filtered = load_filtered()

    df_all = build_dataframe(all_candidates)
    df_top10 = build_filtered_df(filtered)

    print(f"Loaded {len(df_all)} total designs, {len(df_top10)} filtered candidates")
    print(f"Types: {df_all['binder_type'].value_counts().to_dict()}")
    print(f"With structures: {df_all['has_structure'].sum()}")

    # Part 1: Metric distributions
    df_pred = analyze_distributions(df_all, df_top10)
    plot_distributions(df_pred, df_top10)

    # Part 2: Sequence diversity
    identity_matrix, id_df = analyze_diversity(df_top10)
    plot_diversity(df_top10, identity_matrix)

    # Part 3: Scaling analysis
    analyze_scaling(df_all)
    plot_scaling(df_all)

    # Save summary JSON
    summary = {
        "total_designs": len(df_all),
        "with_structures": int(df_all["has_structure"].sum()),
        "top10_count": len(df_top10),
        "type_breakdown": df_all["binder_type"].value_counts().to_dict(),
        "metrics": {},
    }

    df_pred = df_all[df_all["has_structure"]]
    for col in ["iptm", "interface_area", "num_contacts", "ptm", "plddt"]:
        vals = df_pred[col].dropna()
        summary["metrics"][col] = {
            "mean": float(vals.mean()),
            "std": float(vals.std()),
            "median": float(vals.median()),
            "min": float(vals.min()),
            "max": float(vals.max()),
            "q25": float(vals.quantile(0.25)),
            "q75": float(vals.quantile(0.75)),
        }

    # Top 10 summary
    summary["top10"] = []
    for _, row in df_top10.iterrows():
        summary["top10"].append({
            "design_id": row["design_id"],
            "binder_type": row["binder_type"],
            "composite_score": float(row["composite_score"]),
            "interface_area": float(row["interface_area"]),
            "num_contacts": int(row["num_contacts"]),
            "humanness": float(row["humanness"]),
        })

    # Diversity summary
    off_diag = []
    n = len(df_top10)
    for i, j in combinations(range(n), 2):
        off_diag.append(identity_matrix[i, j])
    summary["diversity"] = {
        "mean_pairwise_identity": float(np.mean(off_diag)),
        "min_pairwise_identity": float(np.min(off_diag)),
        "max_pairwise_identity": float(np.max(off_diag)),
    }

    summary_path = OUT_DIR / "analysis_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved: {summary_path}")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nOutputs in: {OUT_DIR}")
    print(f"  - metric_distributions.png")
    print(f"  - sequence_diversity.png")
    print(f"  - scaling_analysis.png")
    print(f"  - analysis_summary.json")


if __name__ == "__main__":
    main()
