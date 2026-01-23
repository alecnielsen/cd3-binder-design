# Handoff Notes - Modal Integration

## Current Status (2026-01-22)

### Boltz-2: WORKING

Successfully tested structure prediction:
```bash
modal run modal/boltz2_app.py --binder-seq "EVQLVES..." --target-seq "DGNE..."
```

Results:
- pTM: 0.770
- ipTM: 0.528
- pLDDT: 0.9 (note: scale seems off, investigate)
- Interface contacts: 0 (mmCIF parsing not implemented, see below)

**Known limitation:** Interface metrics return 0 because mmCIF parsing isn't implemented. The confidence scores from Boltz-2's JSON output work fine.

### BoltzGen: NEEDS INVESTIGATION

The pipeline runs through design and inverse_folding steps but:
1. **Inverse folding doesn't assign sequences** - Output CIF shows chain B still has all X's
2. **Folding step crashes** - TypeError in confidence module

#### What We Observed

Design spec being generated:
```yaml
entities:
  - protein:
      id: A
      sequence: QTPYKVSISGTTVILTCPQYPGSEILWQHNDKNIGGDEDDKNIGSDEDHLSLKEFSELEQSGYYVCYPRGSKPEDANFYLYLRARVCENCM
      design: false
  - protein:
      id: B
      sequence: XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      design: true
      binding: [A]
```

Output CIF after "successful" inverse_folding step:
```
1 polypeptide(L) ? QTPYKVSISGTTVILTCPQYPGSEILWQHNDKNIGGDEDDKNIGSDEDHLSLKEFSELEQS  <- Target OK
2 polypeptide(L) ? XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX   <- Binder still X's!
```

#### Likely Issues to Investigate

1. **Design spec format may be wrong** - Our YAML structure (`design: true`, `binding: [A]`) may not match BoltzGen's expected input format. Check:
   - BoltzGen GitHub examples
   - OpenProtein's implementation (user says they have it working)
   - `boltzgen run --help` for input format details

2. **Protocol selection** - We're using `--protocol nanobody-anything`. Other options exist:
   - `protein-anything`
   - `peptide-anything`
   - `antibody-anything`

3. **Input might need PDB structure, not just sequence** - The design step might need the target's 3D structure, not just sequence in YAML

4. **Binder specification** - The `sequence: XXX...` placeholder might not be the right way to specify "design this chain"

#### Files Changed

- `modal/boltzgen_app.py`: Added `from __future__ import annotations` for Python 3.9 compatibility
- `modal/boltzgen_app.py`: Fixed download to use `boltzgen/boltzgen-1` repo (all checkpoints in one repo)
- `modal/boltzgen_app.py`: Fixed inference-data download to use `repo_type="dataset"`

### Model Downloads: COMPLETE

Both model weight downloads work:
```bash
modal run modal/boltz2_app.py --download      # ~2GB from boltz-community/boltz-2
modal run modal/boltzgen_app.py --download    # ~10GB from boltzgen/boltzgen-1 + boltzgen/inference-data
```

Models stored in Modal volumes:
- `boltz-models` - Boltz-2 weights
- `boltzgen-models` - BoltzGen weights

## Next Steps for Investigation

1. **Check BoltzGen input format**
   - Look at examples in https://github.com/HannesStark/boltzgen
   - The input might need a PDB file for the target structure, not just sequence
   - Check if there's a specific format for specifying designable regions

2. **Compare with OpenProtein implementation**
   - User mentioned OpenProtein has BoltzGen working
   - Find their implementation and compare input formats

3. **Test with official BoltzGen examples**
   - Clone the BoltzGen repo
   - Run their example notebooks/scripts
   - Compare our input format with theirs

4. **Check if target needs structure**
   - BoltzGen might need the target PDB structure for spatial reasoning
   - Try passing actual PDB file instead of extracting sequence

## Technical Notes

### Modal CLI Location
```bash
/Users/alec/Library/Python/3.9/bin/modal
```

### Test Commands
```bash
# Boltz-2 (working)
/Users/alec/Library/Python/3.9/bin/modal run modal/boltz2_app.py --binder-seq "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDYWMNWVRQAPGKGLEWVAEIRLKSNNYATHYAESVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCTGSYYGMDYWGQGTLVTVSS" --target-seq "DGNEEMGGITQTPYKVSISGTTVILTCPQYPGSEILWQHNDKNIGGDEDDKNIGSDEDHLSLKEFSELEQSGYYVCYPRGSKPEDANFYLYLRARVCENCME"

# BoltzGen (not working correctly)
/Users/alec/Library/Python/3.9/bin/modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb --target-chain A --num-designs 2
```

### Key Error in BoltzGen Folding Step
```
TypeError: expected Tensor as element 0 in argument 0, but got float
```
At: `/usr/local/lib/python3.11/site-packages/boltzgen/model/modules/confidence.py:143`

This might be a secondary issue - the primary problem is that inverse folding isn't assigning sequences.

## References

- BoltzGen GitHub: https://github.com/HannesStark/boltzgen
- BoltzGen HuggingFace: https://huggingface.co/boltzgen/boltzgen-1
- Boltz GitHub: https://github.com/jwohlwend/boltz
- Working Boltz implementation: `/Users/alec/kernel/protein-predict/boltz_modal.py`
