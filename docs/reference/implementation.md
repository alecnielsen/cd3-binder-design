# Implementation Details

Code examples and implementation patterns used in the pipeline.

## Configuration

Pipeline behavior is controlled by `config.yaml`:

```yaml
design:
  num_vhh_designs: 200          # VHH nanobodies via nanobody-anything
  num_fab_designs: 200          # Fab CDR redesign via antibody-anything
  target_structures:
    - data/targets/1XIW.pdb     # CD3εδ heterodimer
    - data/targets/1SY6.pdb     # CD3εγ heterodimer

  # Human Fab scaffolds for CDR redesign
  fab_scaffolds:
    - adalimumab                # 6cr1
    - belimumab                 # 5y9k
    - dupilumab                 # 6wgb
  fab_scaffold_dir: data/fab_scaffolds

  # Optimization track (existing binders)
  starting_sequences:
    - teplizumab
    - sp34
    - ucht1
  affinity_variants:
    - wild_type
    - 10x_weaker
    - 100x_weaker

calibration:
  positive_controls:
    - teplizumab
    - sp34
    - ucht1
  calibration_margin:
    pdockq: 0.05
    interface_area: 100
    contacts: 2

filtering:
  binding:
    min_pdockq: 0.5
    min_interface_area: 800
    min_contacts: 10
    use_calibrated: true

  humanness:
    min_oasis_score: 0.8
    generate_back_mutations: true

  liabilities:
    allow_deamidation_cdr: false
    allow_isomerization_cdr: false
    allow_glycosylation_cdr: false
    max_oxidation_sites: 2

  developability:
    cdr_h3_length_range: [8, 20]
    net_charge_range: [-2, 4]
    pi_range: [6.0, 9.0]
    max_hydrophobic_patches: 2

  fallback:
    min_candidates: 10
    relax_soft_filters_first: true
    max_threshold_relaxation: 0.1

formatting:
  tumor_target: trastuzumab
  formats:
    - crossmab
    - fab_scfv
    - fab_vhh
    - igg_scfv
    - igg_vhh
  linkers:
    scfv: "GGGGSGGGGSGGGGS"
    fc_fusion: "GGGGSGGGGS"

output:
  num_final_candidates: 10
  include_structures: true
  generate_report: true

reproducibility:
  boltzgen_seed: 42
  sampling_seed: 12345
  clustering_seed: 0
```

## Sequence Liability Detection

```python
# src/analysis/liabilities.py

DEAMIDATION_MOTIFS = ['NG', 'NS', 'NT', 'ND', 'NH']
ISOMERIZATION_MOTIFS = ['DG', 'DS', 'DT', 'DD', 'DH', 'DN']
OXIDATION_RESIDUES = ['M', 'W']

def find_glycosylation_sites(sequence: str) -> list[int]:
    """Find N-X-S/T motifs where X != P."""
    sites = []
    for i in range(len(sequence) - 2):
        if (sequence[i] == 'N' and
            sequence[i+1] != 'P' and
            sequence[i+2] in ['S', 'T']):
            sites.append(i)
    return sites

def count_unpaired_cysteines(sequence: str) -> int:
    """Odd cysteine count suggests unpaired."""
    return sequence.count('C') % 2
```

## Humanness Scoring

```python
# src/analysis/humanness.py
from biophi import Sapiens, OASis

def score_humanness(vh_sequence: str, vl_sequence: str = None) -> dict:
    """Score humanness using BioPhi OASis."""
    oasis = OASis()
    vh_score = oasis.score(vh_sequence, chain_type='H')
    result = {'vh_oasis': vh_score}

    if vl_sequence:
        vl_score = oasis.score(vl_sequence, chain_type='L')
        result['vl_oasis'] = vl_score
        result['mean_oasis'] = (vh_score + vl_score) / 2

    return result
```

## Bispecific Assembly

```python
# src/formatting/crossmab.py

def assemble_crossmab(
    target_vh: str, target_vl: str,
    cd3_vh: str, cd3_vl: str,
    fc_knob: str, fc_hole: str,
    ch1: str, cl: str, hinge: str
) -> dict:
    """Assemble CrossMab bispecific sequences."""
    # Target arm: standard Fab
    heavy_chain_1 = target_vh + ch1 + hinge + fc_knob
    light_chain_1 = target_vl + cl

    # CD3 arm: CrossMab swap (CH1-CL exchanged)
    heavy_chain_2 = cd3_vh + cl + hinge + fc_hole  # VH fused to CL
    light_chain_2 = cd3_vl + ch1                    # VL fused to CH1

    return {
        'heavy_chain_1': heavy_chain_1,
        'heavy_chain_2': heavy_chain_2,
        'light_chain_1': light_chain_1,
        'light_chain_2': light_chain_2,
    }
```

## Understanding CIF Chain IDs

**Critical for Fab scaffolds**: BoltzGen uses CIF chain identifiers (`struct_asym.id`) which differ from PDB author chain IDs.

### CIF vs PDB Chain Identifiers

| Identifier | Source | Example | Meaning |
|------------|--------|---------|---------|
| `struct_asym.id` | CIF format, BoltzGen | A, B, C, D | Sequential IDs for each polymer entity |
| `pdbx_strand_id` | RCSB API, PDB format | H, L | Author-assigned labels (H=heavy, L=light) |

**Key point**: A, B, C, D are just sequential identifiers assigned to each entity in order of appearance. They do NOT inherently mean "heavy" or "light" chain.

### How CIF Chain IDs Are Assigned

```
PDB Structure Deposited:
  Entity 1: Light chain      → struct_asym.id = A
  Entity 2: Heavy chain      → struct_asym.id = B
  Entity 3: Antigen (target) → struct_asym.id = C
  Entity 4: Ligand           → struct_asym.id = D
```

### Example: Adalimumab (6CR1)

```
Entity 1: Light chain (VL domain) → Chain A in CIF
Entity 2: Heavy chain (VH domain) → Chain B in CIF
```

So the YAML scaffold uses: `vh_chain: "B"`, `vl_chain: "A"`

### Verifying Chain IDs in CIF Files

```bash
# Check struct_asym section (maps chain ID to entity)
grep "_struct_asym" scaffold.cif -A 10 | grep -E "^[A-Z]"
# Output: A N N 1 ?   <- Chain A = Entity 1
#         B N N 2 ?   <- Chain B = Entity 2

# Check entity descriptions
grep -A 200 "^_entity_poly.entity_id" scaffold.cif | head -30
# Look for "heavy chain" or "light chain" in descriptions
```

### All Verified Scaffold Chain Mappings

| Scaffold | PDB | VH | VL | Notes |
|----------|-----|----|----|-------|
| adalimumab | 6cr1 | B | A | VL=entity1, VH=entity2 |
| belimumab | 5y9k | B | A | VL=entity1, VH=entity2 |
| crenezumab | 5vzy | A | B | VH=entity1, VL=entity2 |
| dupilumab | 6wgb | A | B | VH=entity1, VL=entity2 |
| golimumab | 5yoy | G | D | TNF=entity1, VL=entity2, VH=entity3 |
| guselkumab | 4m6m | B | A | VL=entity1, VH=entity2 |
| mab1 | 3h42 | D | C | target=entities1-2, VL=entity3, VH=entity4 |
| necitumumab | 6b3s | B | C | EGFR=entity1, VH=entity2, VL=entity3 |
| nirsevimab | 5udc | A | B | VH=entity1, VL=entity2 |
| sarilumab | 8iow | D | C | IL6R=entity1, VL=entity2, VH=entity3 |
| secukinumab | 6wio | A | B | VH=entity1, VL=entity2 |
| tezepelumab | 5j13 | C | B | TSLP=entity1, VL=entity2, VH=entity3 |
| tralokinumab | 5l6y | B | C | IL13=entity1, VH=entity2, VL=entity3 |
| ustekinumab | 3hmw | B | A | VL=entity1, VH=entity2 |

---

## BoltzGen Fab Scaffold Format

Fab CDR redesign uses YAML scaffold specifications that define which regions to redesign:

```yaml
# data/fab_scaffolds/adalimumab.6cr1.yaml
path: adalimumab.6cr1.cif

include:
  - chain:
      id: B                     # VH chain (entity 2 in this CIF)
      res_index: 1..121
  - chain:
      id: A                     # VL chain (entity 1 in this CIF)
      res_index: 1..107

design:                         # CDR regions to redesign
  - chain:
      id: B
      res_index: 26..32,52..57,99..110  # H-CDR1, H-CDR2, H-CDR3
  - chain:
      id: A
      res_index: 24..34,50..56,89..97   # L-CDR1, L-CDR2, L-CDR3

design_insertions:              # Variable CDR lengths
  - insertion:
      id: B
      res_index: 99
      num_residues: 3..21       # H-CDR3 can vary 3-21 aa
  - insertion:
      id: A
      res_index: 89
      num_residues: 3..12       # L-CDR3 can vary 3-12 aa
```

Setup scaffolds with:
```bash
python scripts/setup_fab_scaffolds.py
```

Verify chain IDs with:
```bash
python scripts/verify_fab_chains.py
```

## Modal Deployment for BoltzGen

```python
# modal/boltzgen_app.py - VHH design
@app.function(gpu="A100", timeout=3600)
def run_boltzgen(
    target_pdb_content: str,
    target_chain: str = "A",
    num_designs: int = 100,
    binder_length: int = 120,
    protocol: str = "nanobody-anything",
    seed: int = 42,
) -> list[dict]:
    """Run BoltzGen VHH design."""
    ...

# modal/boltzgen_app.py - Fab CDR redesign
@app.function(gpu="A100", timeout=5400)
def run_boltzgen_fab(
    target_pdb_content: str,
    target_chain: str = "A",
    scaffold_names: list[str],        # ["adalimumab", "belimumab"]
    scaffold_yaml_contents: dict,      # {name: yaml_content}
    scaffold_cif_contents: dict,       # {name: cif_content}
    num_designs: int = 100,
    seed: int = 42,
) -> list[dict]:
    """Run BoltzGen Fab CDR redesign with human scaffolds."""
    # Returns list with vh_sequence, vl_sequence, pTM, ipTM
    ...
```

## Expected Output Format

Per-candidate JSON scorecard:

```json
{
  "candidate_id": "denovo_vhh_042",
  "source": "boltzgen_1XIW",
  "binding": {
    "iptm": 0.72,
    "iptm_note": "Interface structural confidence, NOT affinity predictor",
    "interface_area": 1245.3,
    "contacts": 18,
    "passed_calibrated_threshold": true
  },
  "epitope": {
    "cd3_contact_residues": [23, 25, 26, 28, 29, 31, 45, 47],
    "okt3_overlap_percent": 0.75,
    "epitope_class": "OKT3-like"
  },
  "humanness": {
    "oasis_score": 0.87,
    "sapiens_mutations_applied": 3,
    "back_mutation_variant_generated": true
  },
  "liabilities": {
    "deamidation_sites": [],
    "glycosylation_sites": [],
    "isomerization_sites": [],
    "oxidation_sites": [{"position": 55, "residue": "M", "severity": "soft"}],
    "unpaired_cys": 0
  },
  "developability": {
    "cdr_h3_length": 14,
    "net_charge": 1,
    "isoelectric_point": 7.8,
    "hydrophobic_patches": 1
  },
  "composite_score": 0.82,
  "rank": 3,
  "_provenance": {
    "pipeline_version": "1.0.0",
    "git_commit": "abc1234",
    "run_timestamp": "2026-01-15T10:30:00Z"
  }
}
```

## Reproducibility

### Version Pinning

All tool versions are pinned in `requirements.txt`:

```
boltzgen==0.2.1
boltz2==0.3.0
anarci==1.3.1
biophi==1.0.5
abnumber==0.3.5
abodybuilder2==1.0.2
```

### Output Provenance

Each output file includes provenance metadata:

```yaml
_provenance:
  pipeline_version: "1.0.0"
  git_commit: "abc1234"
  run_timestamp: "2026-01-15T10:30:00Z"
  config_hash: "sha256:..."
  tool_versions:
    boltzgen: "0.2.1"
    boltz2: "0.3.0"
```

## Computational Resources

| Step | GPU Type | Time (est.) | Cost (est.) |
|------|----------|-------------|-------------|
| BoltzGen (400 designs) | A100 | 1-2 hours | $5-15 |
| Boltz-2 (400 complexes) | A100 | 2-4 hours | $10-30 |
| ABodyBuilder2 | A100 | 30 min | $2-5 |
| **Total** | | **4-8 hours** | **$20-50** |
