# Phase 3 Enhancement Summary: Key Insights & Quick Reference

## Problem Statement

Current implementation (Phase 3 initial) correctly applies the SSH step's mix numbering strategy to the assembly model, but **TDT enforcement only happens for the SSH step**. Flux steps reuse SSH mix numbering even if they specify a different strategy in the YAML.

**Impact**: Users cannot optimize flux step tracking with a different mix granularity (e.g., finer by_pin tracking in L1, coarser by_material in L2).

---

## Solution Architecture

### Three-Tier Design

```
┌─────────────────────────────────────────────────────────────┐
│ 1. ASSEMBLY LAYER (DragonModel.py)                          │
│    • Multi-strategy state machine                           │
│    • Strategy switching with correspondence tracking        │
│    • Per-step mix state snapshots                           │
└─────────────────────────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────┐
│ 2. CASE GENERATION LAYER (case_generator.py)               │
│    • Per-step strategy application (refactor needed)       │
│    • Per-step TDT enforcement (new)                        │
│    • Correspondence recording between steps                │
└─────────────────────────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────┐
│ 3. PROCEDURE LAYER (dragon_module_calls.py)                │
│    • TRK.c2m: per-step mix references (enhance)           │
│    • EDIR.c2m: per-step mix references (enhance)          │
│    • Library (MIX.c2m): SSH mixes only (invariant)        │
└─────────────────────────────────────────────────────────────┘
```

### Critical Design Principle: Library Invariance

```
┌──────────────────────────────────────────────────────────────┐
│  MIX.c2m (DRAGON5 LIB: module)                               │
│  ├─ ALWAYS uses SSH step mix numbering                       │
│  ├─ Contains isotopic definitions for ALL fuel mixes        │
│  ├─ Cannot change during execution                          │
│  └─ Acts as the "ground truth" for mix definitions          │
└──────────────────────────────────────────────────────────────┘
       ↑ SSH defines library mixes (immutable)
       │
       ├─→ TRK.c2m (SSH tracking)
       ├─→ TRK.c2m (L1 tracking) ← Can use different mixes
       ├─→ TRK.c2m (L2 tracking) ← Can use different mixes
       │
       └─ EDIR.c2m (outputs) ← Can reference per-step mixes
```

**Key Insight**: Library must be defined once with SSH mixes. Flux steps that use different mix numbering strategies must maintain correspondence mappings back to SSH mixes for material property lookups and depletion tracking.

---

## Workflow Comparison: Before vs After

### Before (Current Phase 3)
```
SSH step:
  1. Apply radii
  2. Apply strategy (by_material)
  3. Enforce TDT

FLUX steps:
  1. Use same mixes as SSH
  (No per-step strategy or TDT enforcement)
```

### After (Enhanced Phase 3)
```
SSH step:
  1. Apply radii
  2. Apply strategy (by_material) → Snapshot + history
  3. Enforce TDT → Map correspondence if needed
  4. Generate library (MIX.c2m) — final indices

FLUX_L1 step:
  1. Check if strategy differs from SSH
  2. If different: apply new strategy (by_pin) → Snapshot + correspondence
  3. Build geometry with glow
  4. Enforce L1-specific TDT indices → Re-map correspondence
  5. Store step_mix_mapping["FLUX_L1"] = {...}

FLUX_L2 step:
  1. Check if strategy differs from SSH/L1
  2. If different: apply new strategy (by_material) → Snapshot + correspondence
  3. Build geometry with glow
  4. Enforce L2-specific TDT indices → Re-map correspondence
  5. Store step_mix_mapping["FLUX_L2"] = {...}

Procedure generation:
  1. Generate TRK.c2m with step_mix_mapping (per-step tracking mixes)
  2. Generate EDIR.c2m with step_mix_mapping (per-step output mixes)
  3. MIX.c2m unchanged (SSH mixes only)
```

---

## Implementation Priority & Complexity

### Phase 3A: Data Flow (MEDIUM complexity)
**Effort**: ~6-8 hours

1. Refactor `generate_cle2000_procedures()` to loop over steps sequentially
2. Extract TDT enforcement into parametrized method
3. Collect per-step mix states into `step_mix_mapping`

**Files**: `case_generator.py`

### Phase 3B: Procedure Layer (MEDIUM complexity)
**Effort**: ~4-6 hours

1. Add `step_mix_overrides` parameter to TRK and EDI_COMPO
2. Update reference logic to look up correct mixes per step
3. Maintain backward compatibility

**Files**: `dragon_module_calls.py`

### Phase 3C: Assembly Model Enhancement (LOW complexity)
**Effort**: ~2-3 hours

1. Add `get_step_specific_mixes(step_id)` query method
2. Ensure correspondence mappings properly indexed by step

**Files**: `DragonModel.py`

---

## Test Strategy: Four Tiers

```
┌─ UNIT (Assembly Model) ──────────────────────────────────────┐
│  ✓ apply_mix_numbering_strategy() correctness               │
│  ✓ _capture_mix_state() snapshot accuracy                   │
│  ✓ record_mix_correspondence_transition() mapping logic     │
│  ✓ enforce_tdt_with_correspondence_remapping() integration  │
└──────────────────────────────────────────────────────────────┘
           ↓
┌─ INTEGRATION (YAML → Assembly) ──────────────────────────────┐
│  ✓ Single strategy (backward compat)                        │
│  ✓ Two-level mixed strategies ⭐ PRIMARY                    │
│  ✓ Multi-level correspondence chains                        │
│  ✓ Per-step mix index queries                               │
└──────────────────────────────────────────────────────────────┘
           ↓
┌─ END-TO-END (Full Workflow) ──────────────────────────────────┐
│  ✓ DragonCase with per-step YAML → procedures generated     │
│  ✓ Library preserved across steps                           │
│  ✓ TRK/EDIR updated with per-step references               │
│  ✓ No regressions in single-strategy cases                  │
└──────────────────────────────────────────────────────────────┘
           ↓
┌─ REFERENCE DATA (User Provides) ──────────────────────────────┐
│  • TDT files for multi-strategy scenarios                    │
│  • Test data: by_material_ssh.dat, by_pin_l1.dat, etc.     │
│  • Used by all integration & end-to-end tests              │
└──────────────────────────────────────────────────────────────┘
```

---

## Code Locations & References

### Core Implementation Files

| File | Purpose | Change Type |
|------|---------|-------------|
| `DragonModel.py` | Assembly state machine | ✅ Done (Phase 2) |
| `case_generator.py` | Workflow orchestration | ⏳ Refactor needed |
| `dragon_module_calls.py` | TRK/EDIR/LIB generation | ⏳ Enhancement needed |
| `DragonCalculationScheme.py` | Step definition | ✅ Done (Phase 1) |

### Test Files

| File | Purpose | Status |
|------|---------|--------|
| `tests/test_assembly_model.py` | Unit tests | ⏳ Add mix numbering tests |
| `tests/test_dragon_scheme_generation_mix_numbering.py` | Integration tests | ⏳ Create new |
| `tests/test_dragon_case_per_step_mix_numbering.py` | End-to-end tests | ⏳ Create new |
| `tests/reference_tdt_files/` | Reference data | ⏳ User provides |

---

## Key Methods & Data Structures

### Assembly Model State Machine

```python
# Current active state
assembly.current_mix_numbering_strategy = "by_material"  # or "by_pin"
assembly.fuel_material_mixture_names = ["UOX_zone_1", "UOX_zone_2", ...]
assembly.fuel_material_mixture_indices = [1, 2, 3, ...]

# History & correspondence
assembly.mix_state_history = [
    {"strategy": "by_material", "mix_dict": {...}, "timestamp": 0},
    {"strategy": "by_pin", "mix_dict": {...}, "timestamp": 1},
    ...
]

assembly.mix_correspondence_mapping = {
    ("by_material", "by_pin"): {1: 1, 2: 2, 3: 9, 4: 10},  # temp → temp
    ("by_pin", "by_material"): {1: 1, 2: 2, 9: 3, 10: 4},  # reverse
}

assembly.temp_to_tdt_index_mapping = {
    1: 101, 2: 102, 3: 103, ...  # temp → final TDT index
}
```

### Step Mix Mapping Data Structure

```python
step_mix_mapping = {
    "SSH": {
        "strategy": "by_material",
        "indices": [1, 2, 3, 4],
        "names": ["UOX_zone_1", "UOX_zone_2", "GD_zone_1", "GD_zone_2"],
        "tdt_enforced": True,
    },
    "FLUX_L1": {
        "strategy": "by_pin",
        "indices": [1, 2, ..., 16],
        "names": ["UOX_zone_1_pin1", ..., "GD_zone_2_pin8"],
        "tdt_enforced": True,
        "correspondence": ("by_material", "by_pin"),
    },
    "FLUX_L2": {
        "strategy": "by_material",
        "indices": [1, 2, 3, 4],
        "names": ["UOX_zone_1", "UOX_zone_2", "GD_zone_1", "GD_zone_2"],
        "tdt_enforced": True,
        "correspondence": ("by_pin", "by_material"),
    },
}

# Used by:
trk = TRK(scheme, case_name, step_mix_overrides=step_mix_mapping)
edi_compo = EDI_COMPO(assembly, step_mix_overrides=step_mix_mapping)
```

---

## Example Test Cases

### T1: Unit Test
```python
def test_record_mix_correspondence_transition():
    # Setup
    assembly = CartesianAssemblyModel(...)
    assembly.set_material_compositions(compositions)
    assembly.analyze_lattice_description(build_pins=True)
    
    # by_material (SSH)
    assembly.apply_mix_numbering_strategy("by_material")
    # mix_state_history[-1] = {strategy: "by_material", mix_dict: {1,2,3,4}, ...}
    
    # by_pin (FLUX_L1)
    assembly.apply_mix_numbering_strategy("by_pin")
    # mix_state_history[-1] = {strategy: "by_pin", mix_dict: {1..16}, ...}
    # mix_correspondence_mapping[("by_material", "by_pin")] = {1:1, 2:2, 3:9, 4:10}
    
    # Verify
    assert len(assembly.mix_state_history) == 2
    assert assembly.mix_correspondence_mapping[("by_material", "by_pin")] == {1:1, 2:2, 3:9, 4:10}
```

### T2: Integration Test
```python
def test_two_level_mixed_strategies():
    # YAML with SSH by_material, L1 by_pin, L2 by_material
    case = DragonCase(
        case_name="AT10_hybrid",
        call_glow=False,
        ...
        tdt_path="tests/reference_tdt_files",
    )
    
    # Generate procedures
    result = case.generate_cle2000_procedures()
    
    # Verify assembly state after generation
    assert case.assembly.current_mix_numbering_strategy == "by_material"  # Back to SSH
    assert len(case.assembly.mix_state_history) == 3  # SSH, L1, L2
    assert ("by_material", "by_pin") in case.assembly.mix_correspondence_mapping
    assert ("by_pin", "by_material") in case.assembly.mix_correspondence_mapping
    
    # Verify procedures reference correct mixes
    # (Check TRK.c2m and EDIR.c2m contain expected mix names/indices)
```

### T3: End-to-End Test
```python
def test_dragon_case_preserves_library_across_steps():
    case = DragonCase(..., config_yamls={
        "CALC_SCHEME": "configs/hybrid_mix_numbering.yaml",
        ...
    })
    
    result = case.generate_cle2000_procedures()
    
    # Verify MIX.c2m uses SSH mixes only
    mix_content = read_file(result["mix"])
    assert "UOX_zone_1" in mix_content
    assert "UOX_zone_1_pin1" not in mix_content  # No L1 pin mixes in library
    
    # Verify TRK.c2m has per-step tracking
    trk_content = read_file(result["trk"])
    # L1 tracking should reference by_pin mixes
    # L2 tracking should reference by_material mixes
```

---

## Risk Analysis & Mitigation

| Risk | Impact | Mitigation |
|------|--------|-----------|
| Library mix invariance violated | CRITICAL | Immutable MIX.c2m generation, test validation |
| Correspondence mapping incomplete | HIGH | Comprehensive history tracking, unit tests |
| Backward compat broken | HIGH | Single-strategy cases must produce identical output |
| TDT file not found for flux step | MEDIUM | Clear error messages, optional enforcement |
| Mix index misalignment in EDIR | MEDIUM | Reference data validation, procedure review |

---

## Success Metrics

- ✅ All existing tests continue to pass (backward compatibility)
- ✅ New unit tests for mix numbering state machine (20+ tests)
- ✅ New integration tests for per-step strategies (15+ tests)
- ✅ Generated procedures match reference outputs (visual inspection)
- ✅ Correspondence mappings correctly track all strategy transitions
- ✅ Documentation updated with YAML examples
- ✅ Performance: No significant slowdown (<5% overhead)

---

## Timeline Estimate

| Phase | Task | Effort | Duration |
|-------|------|--------|----------|
| 3A | Data flow refactoring | 6-8h | 1-2 days |
| 3B | Procedure layer enhancement | 4-6h | 1 day |
| 3C | Assembly model extension | 2-3h | 0.5 day |
| Test | Unit + Integration + E2E | 8-10h | 1-2 days |
| **Total** | | **20-27h** | **3-4 days** |

*(Assumes reference TDT data provided by user)*

