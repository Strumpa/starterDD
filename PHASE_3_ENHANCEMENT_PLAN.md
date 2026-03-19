# Phase 3 Enhancement Plan: Per-Step Mix Numbering with Per-Step TDT Enforcement

## Executive Summary

Current Phase 3 implementation (hybrid mix numbering orchestration) correctly applies the SSH step's mix numbering strategy during assembly model initialization, but only enforces TDT indices from the SSH step's TDT file. The enhancement enables **per-step mix numbering strategies** where each calculation step (SSH, FLUX_L1, FLUX_L2) can have its own mix numbering strategy and corresponding TDT file for mix index enforcement.

---

## Current Workflow Analysis

### Current Flow (Lines 245-320 in case_generator.py)

```
1. Load scheme from YAML
2. Build assembly model (lattice, pins, material compositions)
3. Apply SSH step radii subdivision
4. Apply SSH step mix numbering strategy → numbers all mixes according to SSH strategy
5. Loop over all trackable steps:
   - Build geometry with glow (if call_glow=True)
   - If step is SSH: read TDT file & enforce indices
   - If step is FLUX: no TDT enforcement (uses SSH indices)
6. Identify generating/daughter mixes (based on final SSH indices)
7. Generate CLE2000 procedures using SSH mix numbering strategy
```

### Problem Identified

- **SSH Step**: Mix numbering strategy applied ✅, TDT indices enforced ✅
- **FLUX Steps**: No per-step mix numbering strategy support ❌, no per-step TDT enforcement ❌
- **Result**: All steps use SSH's mix numbering and indices, losing per-step optimization opportunity

---

## Proposed Enhancement: Per-Step Mix Numbering Workflow

### New Flow Architecture

```
1. Load scheme from YAML
2. Build assembly model
3. Apply SSH step radii subdivision
4. Apply SSH step mix numbering strategy → numbers all mixes
5. Enforce SSH TDT indices (read SSH step's TDT file)
6. Identify generating/daughter mixes (SSH mixes)
7. Generate MIX.c2m with SSH mixes ← Library definition (immutable across steps)
8. Loop over flux steps (optional per-step mix numbering):
   a. Check if flux step has a different mix_numbering_strategy
   b. If yes AND step has dedicated TDT file:
      - Save current assembly mix state (assembly.mix_state_history)
      - Apply new strategy: assembly.apply_mix_numbering_strategy(step.strategy)
      - Build geometry (if call_glow)
      - Enforce step-specific TDT indices
      - Record correspondence transition for depletion tracking
      - Update step tracking structures
   c. If no:
      - Use SSH mixes throughout
9. Generate tracking (TRK.c2m) - references per-step mix numbering
10. Generate EDIR (EDI/COMPO) - references per-step mix numbering
```

### Key Design Decisions

1. **Library (MIX.c2m) remains SSH-based**: The isotopic library is defined once using SSH mix numbering. This is DRAGON5 requirement (all fuel mixes must exist in LIB).

2. **Tracking (TRK.c2m) becomes per-step aware**: Different steps can have different mix numbering if their tracking uses different mix granularity.

3. **Correspondence tracking**: Assembly model maintains correspondence mappings between SSH and flux step mixes for depletion continuity tracking.

4. **Backward compatibility**: If all steps use same mix numbering strategy, behavior is identical to current implementation.

---

## Implementation Plan

### Phase 3A: Data Flow Architecture (YAML → Assembly State)

**Objective**: Enable reading per-step mix numbering strategies from YAML and applying them sequentially.

#### 1.1 YAML Enhancement (Already Done ✅)
- ✅ `CalculationStep` has `mix_numbering_strategy` attribute
- ✅ Default: `"by_material"` (most conservative)
- ✅ Can be overridden per-step in YAML

**Example YAML**:
```yaml
DRAGON_CALCULATION_SCHEME:
  name: "hybrid_mix_numbering"
  steps:
    - name: "SSH"
      step_type: "self_shielding"
      mix_numbering_strategy: "by_material"  # Library uses this
      
    - name: "FLUX_L1"
      step_type: "flux"
      flux_level: 1
      mix_numbering_strategy: "by_pin"       # L1 tracking uses by_pin
      
    - name: "FLUX_L2"
      step_type: "flux"
      flux_level: 2
      mix_numbering_strategy: "by_material"  # L2 tracking reuses by_material
```

#### 1.2 Case Generator Enhancement
**File**: `starterDD/InterfaceToDD/case_generator.py`

**Changes**:
1. **Refactor `generate_cle2000_procedures()`** (lines 245-370):
   - Separate SSH initialization from flux step processing
   - Extract per-step TDT enforcement logic
   - Add correspondence recording between steps

2. **Add method `_apply_step_mix_numbering(step, assembly, glow_enabled)`**:
   - Check if step mix_numbering_strategy differs from current
   - If yes: call `assembly.apply_mix_numbering_strategy(step.strategy)`
   - Handle correspondence recording
   - Return updated assembly state

3. **Add method `_enforce_step_tdt_indices(step, assembly, tdt_path, tdt_base_name, call_glow)`**:
   - Identical to current SSH logic but parametrized
   - Read TDT file for any trackable step (SSH, FLUX_L1, FLUX_L2)
   - Call `assembly.enforce_material_mixture_indices_from_tdt(tdt_indices)`
   - Record which step the enforcement applies to

#### 1.3 Assembly Model Enhancement (Already Done ✅)
**File**: `starterDD/DDModel/DragonModel.py`

- ✅ `_capture_mix_state()` — snapshots current state
- ✅ `apply_mix_numbering_strategy(strategy)` — switches strategy + captures
- ✅ `record_mix_correspondence_transition()` — tracks transitions
- ✅ Enhanced `enforce_material_mixture_indices_from_tdt()` with correspondence re-mapping

**New enhancement needed**:
4. **Add method `get_step_specific_mixes(step_name_or_index)`**:
   - Returns mix indices/names for a specific step
   - Looks up in correspondence mappings if step used different strategy
   - Used by TRK/EDIR generation to reference correct mixes

---

### Phase 3B: Tracking & Edition Layer Updates

**Objective**: Make TRK.c2m and EDIR.c2m aware of per-step mix numbering.

#### 2.1 TRK Generation Enhancement
**File**: `starterDD/InterfaceToDD/dragon_module_calls.py` (TRK class)

**Changes**:
1. Add optional parameter `step_mix_overrides: dict[str, list[int]]` to `__init__`
2. Store mapping: `flux_step_name → mix_indices` for that step
3. When writing TRK.c2m, reference correct mixes for each step:
   - SSH tracking → uses SSH mixes
   - FLUX_L1 tracking → uses FLUX_L1 mixes (if different)
   - FLUX_L2 tracking → uses FLUX_L2 mixes (if different)

#### 2.2 EDIR/COMPO Enhancement
**File**: `starterDD/InterfaceToDD/dragon_module_calls.py` (EDI_COMPO class)

**Changes**:
1. Add optional parameter `step_mix_overrides: dict[str, list[int]]`
2. Store per-step mix references
3. When writing EDIR procedures, reference correct mixes for each step

---

### Phase 3C: Integration in case_generator.py

**Refactored `generate_cle2000_procedures()` pseudocode**:

```python
def generate_cle2000_procedures(self):
    # 1. Load & initialize
    scheme = DragonCalculationScheme.from_yaml(self.calc_scheme_yaml)
    assembly, compositions = self._build_assembly_model()
    
    # 2. SSH step: initialize and generate library
    ssh_step = scheme.get_self_shielding_steps()[0]
    ssh_step.apply_radii(assembly)
    assembly.apply_mix_numbering_strategy(ssh_step.mix_numbering_strategy)
    ssh_step_mixes = {
        "indices": assembly.fuel_material_mixture_indices,
        "names": assembly.fuel_material_mixture_names,
    }
    
    # 3. SSH TDT enforcement & library generation
    self._enforce_step_tdt_indices(ssh_step, assembly, ...)
    assembly.identify_generating_and_daughter_mixes()
    lib = LIB(assembly, ...)  # Uses SSH mixes
    
    # 4. Per-step mix numbering (flux steps)
    step_mix_mapping = {}  # step_name → mix indices/names
    for flux_step in scheme.get_flux_steps():
        # 4a. Apply flux step mix numbering (if different)
        if flux_step.mix_numbering_strategy != ssh_step.mix_numbering_strategy:
            self._apply_step_mix_numbering(flux_step, assembly)
            # Record correspondence SSH → FLUX
        
        # 4b. Enforce flux step TDT indices (if TDT file exists)
        if has_tdt_for_step(flux_step):
            self._enforce_step_tdt_indices(flux_step, assembly, ...)
            
        # 4c. Store mix state for this step
        step_mix_mapping[flux_step.name] = {
            "indices": assembly.fuel_material_mixture_indices.copy(),
            "names": assembly.fuel_material_mixture_names.copy(),
            "strategy": flux_step.mix_numbering_strategy,
        }
    
    # 5. Generate tracking & edition with per-step mixes
    trk = TRK(scheme, self.case_name, step_mix_overrides=step_mix_mapping)
    edi_compo = EDI_COMPO(assembly, step_mix_overrides=step_mix_mapping)
    
    # 6. Generate main x2m
    x2m_path = self._build_main_x2m(...)
    
    return {...}
```

---

## Testing Strategy

### Test Categories

#### T1: Unit Tests (Assembly Model Mix Numbering)
**File**: `tests/test_assembly_model.py` (new tests)

1. **test_apply_mix_numbering_strategy_by_material_to_by_pin**
   - Start with `by_material` numbering
   - Switch to `by_pin` numbering
   - Verify: mix names include `_pinN` suffix
   - Verify: correspondence mapping created

2. **test_apply_mix_numbering_strategy_by_pin_to_by_material**
   - Start with `by_pin` numbering
   - Switch to `by_material` numbering
   - Verify: mix names revert to `_zone_N` format
   - Verify: correspondence mapping created

3. **test_capture_mix_state_snapshot**
   - Number mixes with `by_material`
   - Capture state
   - Verify: snapshot contains strategy, mix_dict, timestamp

4. **test_record_mix_correspondence_transition**
   - Apply `by_material` → capture
   - Apply `by_pin` → capture + record correspondence
   - Verify: correspondence maps old temp indices to new temp indices

5. **test_enforce_tdt_with_correspondence_remapping**
   - Number mixes `by_pin` (temporary indices)
   - Record correspondences (hypothetical old → new)
   - Enforce TDT indices
   - Verify: correspondence mappings updated with final TDT indices

#### T2: Integration Tests (YAML → Assembly State)
**File**: `tests/test_dragon_scheme_generation_mix_numbering.py` (new file)

1. **test_single_step_by_material**
   - YAML: SSH with `mix_numbering_strategy: "by_material"`
   - Generate procedures
   - Verify: assembly uses by_material mixes
   - Verify: MIX.c2m contains by_material mix names

2. **test_single_step_by_pin**
   - YAML: SSH with `mix_numbering_strategy: "by_pin"`
   - Generate procedures
   - Verify: assembly uses by_pin mixes
   - Verify: MIX.c2m contains by_pin mix names (with _pinN suffix)

3. **test_two_level_same_strategy**
   - YAML: SSH by_material, FLUX_L1 by_material, FLUX_L2 by_material
   - Generate procedures
   - Verify: all steps use same mixes
   - Verify: correspondence mappings empty/unused

4. **test_two_level_mixed_strategies** ⭐ **PRIMARY HYBRID TEST**
   - YAML: SSH by_material, FLUX_L1 by_pin, FLUX_L2 by_material
   - TDT files provided for all three steps (user provides via reference_tdt_files/)
   - Generate procedures
   - Verify: 
     - MIX.c2m uses SSH by_material mixes (library unchanged)
     - TRK.c2m references FLUX_L1 by_pin mixes in L1 tracking
     - TRK.c2m references FLUX_L2 by_material mixes in L2 tracking
     - Correspondence mappings: SSH ↔ L1 and L1 ↔ L2

5. **test_multi_level_correspondence_chain**
   - YAML: SSH by_material → FLUX_L1 by_pin → FLUX_L2 by_pin
   - Verify: complete correspondence chain maintained

6. **test_mix_indices_representation_at_step**
   - After applying different strategies at different steps
   - Call hypothetical method `assembly.get_step_specific_mixes(step_index)`
   - Verify: returns correct mix indices for that step

#### T3: End-to-End Workflow Tests
**File**: `tests/test_dragon_case_per_step_mix_numbering.py` (new file)

1. **test_dragon_case_with_per_step_yaml**
   - Use YAML with per-step mix numbering
   - Call `DragonCase.generate_cle2000_procedures()`
   - Verify: no errors, all procedures generated
   - Verify: assembly.mix_state_history has multiple entries

2. **test_dragon_case_preserves_library_across_steps**
   - Generate case with mixed strategies
   - Verify: MIX.c2m uses SSH strategy (library invariant)
   - Verify: TRK.c2m and EDIR.c2m updated for per-step strategies

#### T4: TDT Reference Data Tests
**Location**: `tests/reference_tdt_files/`

**Your responsibility** (as noted):
- Generate TDT files for multi-strategy scenarios:
  - `by_material_ssh.dat` — SSH step, by_material strategy
  - `by_pin_l1.dat` — FLUX_L1 step, by_pin strategy
  - `by_material_l2.dat` — FLUX_L2 step, by_material strategy
- Store alongside existing TDT files

**Test usage**:
```python
def test_two_level_mixed_strategies():
    # Uses: reference_tdt_files/by_material_ssh.dat
    #       reference_tdt_files/by_pin_l1.dat
    #       reference_tdt_files/by_material_l2.dat
```

---

## Implementation Checklist

### Phase 3A: Data Flow
- [ ] `case_generator.py`: Refactor `generate_cle2000_procedures()` to handle per-step mix numbering
- [ ] `case_generator.py`: Add `_apply_step_mix_numbering(step, assembly, ...)`
- [ ] `case_generator.py`: Add `_enforce_step_tdt_indices(step, assembly, ...)`
- [ ] `DragonModel.py`: Add `get_step_specific_mixes(step_identifier)` method

### Phase 3B: Tracking & Edition Layer
- [ ] `dragon_module_calls.py` (TRK class): Add `step_mix_overrides` parameter
- [ ] `dragon_module_calls.py` (EDI_COMPO class): Add `step_mix_overrides` parameter
- [ ] Update TRK.write_to_c2m() to reference correct mixes per step
- [ ] Update EDIR.write_to_c2m() to reference correct mixes per step

### Phase 3C: Integration
- [ ] Update `_build_main_x2m()` variants to use step_mix_mapping
- [ ] Verify backward compatibility (all steps same strategy = current behavior)

### Test Implementation
- [ ] Unit tests in `test_assembly_model.py`
- [ ] Integration tests in `test_dragon_scheme_generation_mix_numbering.py`
- [ ] End-to-end tests in `test_dragon_case_per_step_mix_numbering.py`
- [ ] User provides TDT reference data files

---

## Example Scenario: AT10 3x3 with Hybrid Mix Numbering

### Workflow
```yaml
DRAGON_CALCULATION_SCHEME:
  name: "AT10_hybrid_mix_numbering_2L"
  steps:
    - name: "SSH"
      step_type: "self_shielding"
      spatial_method: "CP"
      mix_numbering_strategy: "by_material"     # 2-3 mixes
      
    - name: "FLUX_L1"
      step_type: "flux"
      flux_level: 1
      spatial_method: "CP"
      mix_numbering_strategy: "by_pin"          # 8 unique pin positions
      
    - name: "Edition_L1_to_L2"
      step_type: "edition_between_levels"
      
    - name: "FLUX_L2"
      step_type: "flux"
      flux_level: 2
      spatial_method: "MOC"
      mix_numbering_strategy: "by_material"     # Back to coarse mixes
```

### Resulting State at Key Points

**After SSH initialization**:
```
assembly.current_mix_numbering_strategy = "by_material"
assembly.fuel_material_mixture_names = ["UOX_zone_1", "UOX_zone_2", "GD_zone_1", "GD_zone_2"]
assembly.fuel_material_mixture_indices = [1, 2, 3, 4]
```

**After FLUX_L1 strategy switch**:
```
assembly.current_mix_numbering_strategy = "by_pin"
assembly.fuel_material_mixture_names = ["UOX_zone_1_pin1", "UOX_zone_2_pin1", ..., "GD_zone_1_pin8"]
assembly.fuel_material_mixture_indices = [1, 2, ..., 16]
assembly.mix_correspondence_mapping = {
    ("by_material", "by_pin"): {1: 1, 2: 2, 3: 9, 4: 10}  # Old temp → New temp
}
```

**After FLUX_L2 strategy revert**:
```
assembly.current_mix_numbering_strategy = "by_material"
assembly.fuel_material_mixture_names = ["UOX_zone_1", "UOX_zone_2", "GD_zone_1", "GD_zone_2"]
assembly.fuel_material_mixture_indices = [1, 2, 3, 4]  # After TDT enforcement
assembly.mix_correspondence_mapping = {
    ("by_material", "by_pin"): {1: 1, 2: 2, 3: 9, 4: 10},
    ("by_pin", "by_material"): {1: 1, 2: 2, 9: 3, 10: 4}  # Reverse mapping
}
```

### Generated Procedures
- **MIX.c2m**: Uses SSH's by_material mixes (invariant library)
- **TRK.c2m**: 
  - SSH tracking → by_material mixes (L1 condensation preparation)
  - L1 tracking → by_pin mixes (detailed pin-wise flux solution)
  - L2 tracking → by_material mixes (final homogenized solution)
- **EDIR.c2m**: Outputs at L2 with by_material mixes

---

## Success Criteria

1. ✅ **Data Flow**: YAML per-step mix numbering strategies correctly parsed and applied
2. ✅ **State Tracking**: Assembly model maintains accurate mix_state_history and correspondence_mapping
3. ✅ **TDT Enforcement**: Each step's TDT indices properly enforced and corresponded
4. ✅ **Procedure Generation**: MIX/TRK/EDIR procedures reference correct mixes per step
5. ✅ **Backward Compatibility**: Single-strategy cases produce identical output to Phase 2
6. ✅ **Test Coverage**: All test categories passing with reference TDT data
7. ✅ **Documentation**: Docstrings updated, usage examples provided

---

## Next Steps

1. **You (User)**:
   - Review and approve this plan
   - Generate reference TDT files for test scenarios
   - Provide any YAML template examples for multi-strategy schemes

2. **Agent (Implementation)**:
   - Implement Phase 3A: Data flow refactoring
   - Implement Phase 3B: Tracking layer updates
   - Implement Phase 3C: Integration and backward compatibility
   - Implement test suite

3. **Validation**:
   - Run all tests with reference TDT data
   - Verify procedures generate correctly
   - Compare with AT10 known-good outputs (if available)

