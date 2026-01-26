# Handoff Notes - Modal Integration

## Current Status (2026-01-25)

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

### BoltzGen: WORKING (Fixed 2026-01-25)

Successfully generating de novo binder designs:
```bash
modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb --target-chain A --num-designs 10
```

Example results:
- design_spec_0: 122 aa, ipTM=0.271, pTM=0.713 (VH-like sequence)
- design_spec_1: 116 aa, ipTM=0.157, pTM=0.741 (VH-like sequence)
- design_spec_2: 110 aa, ipTM=0.171, pTM=0.722 (VL-like sequence)

**Note:** BoltzGen's internal filtering (`filter_rmsd < 2.5`) often rejects all designs. We extract designs anyway for our own filtering pipeline.

#### Root Cause (Fixed)

The YAML format was completely wrong. BoltzGen uses a different specification:

**OLD (Wrong):**
```yaml
entities:
  - protein:
      id: A
      sequence: QTPYKVSISGTTVILT...  # Inline sequence
      design: false
  - protein:
      id: B
      sequence: XXXX...  # Placeholder
      design: true
      binding: [A]
```

**NEW (Correct):**
```yaml
entities:
  - protein:
      id: B
      sequence: 110..130  # Range notation for design length
  - file:
      path: target.pdb  # File reference, not inline sequence
      include:
        - chain:
            id: A
```

Key fixes:
1. Use `sequence: 110..130` (range notation) instead of `sequence: XXX...`
2. Reference target via `file: { path: ... }` instead of inline sequence
3. Remove invalid keys (`design: true`, `binding: [A]`)
4. Binder entity comes FIRST, target comes SECOND

#### Files Changed

- `modal/boltzgen_app.py`: Complete rewrite of `build_design_spec_yaml()` function
- `modal/boltzgen_app.py`: Write target PDB to file and reference it in YAML
- `modal/boltzgen_app.py`: Fixed CIF parsing to extract from `_entity_poly.pdbx_seq_one_letter_code`
- `modal/boltzgen_app.py`: Added CSV metrics parsing for ipTM, pTM, etc.

### Model Downloads: COMPLETE

Both model weight downloads work:
```bash
modal run modal/boltz2_app.py --download      # ~2GB from boltz-community/boltz-2
modal run modal/boltzgen_app.py --download    # ~10GB from boltzgen/boltzgen-1 + boltzgen/inference-data
```

Models stored in Modal volumes:
- `boltz-models` - Boltz-2 weights
- `boltzgen-models` - BoltzGen weights

## Next Steps

1. **Run full pipeline** - Both BoltzGen and Boltz-2 are working. Run the full pipeline:
   ```bash
   python scripts/00_run_calibration.py
   python scripts/run_full_pipeline.py --config config.yaml
   ```

2. **Generate more designs** - Current test used 3 designs. Production run should use 100+ per target.

3. **Test with hotspot residues** - Can specify binding site residues for more targeted design:
   ```yaml
   binding_types:
     - chain:
         id: A
         binding: 23,24,25,50,51  # OKT3 epitope residues
   ```

4. **Consider protocol options** - Current uses `nanobody-anything`. Also available:
   - `protein-anything` - general protein binder design
   - `antibody-anything` - for Fab/scFv with CDR design

## Technical Notes

### Modal CLI Location
```bash
/Users/alec/Library/Python/3.9/bin/modal
```

### Test Commands
```bash
# Boltz-2 (working)
modal run modal/boltz2_app.py --binder-seq "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDYWMNWVRQAPGKGLEWVAEIRLKSNNYATHYAESVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCTGSYYGMDYWGQGTLVTVSS" --target-seq "DGNEEMGGITQTPYKVSISGTTVILTCPQYPGSEILWQHNDKNIGGDEDDKNIGSDEDHLSLKEFSELEQSGYYVCYPRGSKPEDANFYLYLRARVCENCME"

# BoltzGen (working)
modal run modal/boltzgen_app.py --target-pdb data/targets/1XIW.pdb --target-chain A --num-designs 10

# Download models (first time only)
modal run modal/boltzgen_app.py --download
modal run modal/boltz2_app.py --download
```

## References

- BoltzGen GitHub: https://github.com/HannesStark/boltzgen
- BoltzGen HuggingFace: https://huggingface.co/boltzgen/boltzgen-1
- Boltz GitHub: https://github.com/jwohlwend/boltz
- Working Boltz implementation: `/Users/alec/kernel/protein-predict/boltz_modal.py`
