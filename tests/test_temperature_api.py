"""
Tests for the new temperature registration API in DragonCase.

Tests cover:
1. API registration (set_fuel_material_temperatures, set_non_fuel_temperatures)
2. Material type key loading from YAML
3. Temperature assignment to MaterialMixture objects
4. Deprecation of associate_temperatures_from_materials_yaml
"""

import pytest
import tempfile
import os
from pathlib import Path
import yaml
import warnings
import sys

from starterDD.MaterialProperties.material_mixture import Composition, MaterialMixture, parse_all_compositions_from_yaml
from starterDD.DDModel import CartesianAssemblyModel
from starterDD.InterfaceToDD.case_generator import DragonCase
from starterDD.DDModel.helpers import associate_temperatures_from_materials_yaml


class TestMaterialTypeKeyLoading:
    """Test loading material_type_key from YAML MATERIAL_TYPES section."""

    def test_material_type_key_loaded_from_yaml(self):
        """Verify material_type_key is set on Composition from YAML MATERIAL_TYPES"""
        yaml_content = {
            'MIX_COMPOSITIONS': [
                {
                    'name': 'UOX_4.5',
                    'isotopic_composition': {'U235': 1e-4, 'U238': 2e-2, 'O16': 4e-2}
                },
                {
                    'name': 'Zircaloy',
                    'isotopic_composition': {'Zr90': 5e-2}
                }
            ],
            'MATERIAL_TYPES': {
                'UOX_4.5': 'fuel',
                'Zircaloy': 'structural'
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(yaml_content, f)
            yaml_file = f.name

        try:
            compositions = parse_all_compositions_from_yaml(yaml_file)

            # Find compositions by name
            uox = next((c for c in compositions if c.material_name == 'UOX_4.5'), None)
            zr = next((c for c in compositions if c.material_name == 'Zircaloy'), None)

            assert uox is not None, "UOX_4.5 composition not found"
            assert zr is not None, "Zircaloy composition not found"

            assert uox.material_type_key == 'fuel', f"Expected 'fuel', got {uox.material_type_key}"
            assert zr.material_type_key == 'structural', f"Expected 'structural', got {zr.material_type_key}"
        finally:
            os.unlink(yaml_file)

    def test_material_type_key_none_if_not_in_yaml(self):
        """Verify material_type_key is None if not in MATERIAL_TYPES"""
        yaml_content = {
            'MIX_COMPOSITIONS': [
                {
                    'name': 'UOX_4.5',
                    'isotopic_composition': {'U235': 1e-4, 'U238': 2e-2, 'O16': 4e-2}
                }
            ],
            'MATERIAL_TYPES': {}  # Empty mapping
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(yaml_content, f)
            yaml_file = f.name

        try:
            compositions = parse_all_compositions_from_yaml(yaml_file)
            uox = compositions[0]
            assert uox.material_type_key is None
        finally:
            os.unlink(yaml_file)


class TestMaterialMixtureTypeKey:
    """Test MaterialMixture material_type_key parameter."""

    def test_material_mixture_stores_type_key(self):
        """Verify MaterialMixture stores material_type_key"""
        comp = Composition('UOX', {'U235': 1e-4})
        mix = MaterialMixture(
            material_name='UOX',
            material_mixture_index=1,
            composition=comp,
            temperature=900.0,
            material_type_key='fuel'
        )

        assert mix.material_type_key == 'fuel'
        assert mix.temperature == 900.0

    def test_material_mixture_type_key_optional(self):
        """Verify material_type_key is optional (defaults to None)"""
        comp = Composition('UOX', {'U235': 1e-4})
        mix = MaterialMixture(
            material_name='UOX',
            material_mixture_index=1,
            composition=comp,
            temperature=900.0
        )

        assert mix.material_type_key is None


class TestDragonCaseTemperatureAPI:
    """Test DragonCase set_*_temperatures() API methods."""

    def test_set_fuel_material_temperatures(self):
        """Verify set_fuel_material_temperatures stores fuel temps"""
        # Create minimal YAML files
        mats_yaml = self._create_minimal_mats_yaml()
        geom_yaml = self._create_minimal_geom_yaml()
        scheme_yaml = self._create_minimal_scheme_yaml()

        try:
            case = DragonCase(
                case_name='test',
                call_glow=False,
                draglib_name_to_alias={'test.lib': 'LIB'},
                config_yamls={
                    'MATS': mats_yaml,
                    'GEOM': geom_yaml,
                    'CALC_SCHEME': scheme_yaml
                }
            )

            # Register fuel temperatures
            case.set_fuel_material_temperatures({
                'UOX_4.5': 900.0,
                'MOX_8.0': 1100.0
            })

            assert case._fuel_material_temperatures['UOX_4.5'] == 900.0
            assert case._fuel_material_temperatures['MOX_8.0'] == 1100.0
        finally:
            self._cleanup_yaml_files(mats_yaml, geom_yaml, scheme_yaml)

    def test_set_non_fuel_temperatures(self):
        """Verify set_non_fuel_temperatures stores non-fuel temps"""
        mats_yaml = self._create_minimal_mats_yaml()
        geom_yaml = self._create_minimal_geom_yaml()
        scheme_yaml = self._create_minimal_scheme_yaml()

        try:
            case = DragonCase(
                case_name='test',
                call_glow=False,
                draglib_name_to_alias={'test.lib': 'LIB'},
                config_yamls={
                    'MATS': mats_yaml,
                    'GEOM': geom_yaml,
                    'CALC_SCHEME': scheme_yaml
                }
            )

            # Register non-fuel temperatures
            case.set_non_fuel_temperatures(
                structural_temperature=625.0,
                gap_temperature=600.0,
                coolant_temperature=560.0,
                moderator_temperature=600.0
            )

            assert case._non_fuel_temperatures['structural'] == 625.0
            assert case._non_fuel_temperatures['gap'] == 600.0
            assert case._non_fuel_temperatures['coolant'] == 560.0
            assert case._non_fuel_temperatures['moderator'] == 600.0
        finally:
            self._cleanup_yaml_files(mats_yaml, geom_yaml, scheme_yaml)

    def test_set_multiple_times_updates(self):
        """Verify calling set_fuel_material_temperatures multiple times updates"""
        mats_yaml = self._create_minimal_mats_yaml()
        geom_yaml = self._create_minimal_geom_yaml()
        scheme_yaml = self._create_minimal_scheme_yaml()

        try:
            case = DragonCase(
                case_name='test',
                call_glow=False,
                draglib_name_to_alias={'test.lib': 'LIB'},
                config_yamls={
                    'MATS': mats_yaml,
                    'GEOM': geom_yaml,
                    'CALC_SCHEME': scheme_yaml
                }
            )

            case.set_fuel_material_temperatures({'UOX_4.5': 900.0})
            case.set_fuel_material_temperatures({'MOX_8.0': 1100.0})

            assert case._fuel_material_temperatures['UOX_4.5'] == 900.0
            assert case._fuel_material_temperatures['MOX_8.0'] == 1100.0
        finally:
            self._cleanup_yaml_files(mats_yaml, geom_yaml, scheme_yaml)

    # Helper methods to create minimal YAML files
    @staticmethod
    def _create_minimal_mats_yaml():
        yaml_content = {
            'MIX_COMPOSITIONS': [
                {
                    'name': 'UOX_4.5',
                    'isotopic_composition': {'U235': 1e-4, 'U238': 2e-2}
                }
            ],
            'MATERIAL_TYPES': {'UOX_4.5': 'fuel'}
        }
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
        yaml.dump(yaml_content, f)
        f.close()
        return f.name

    @staticmethod
    def _create_minimal_geom_yaml():
        yaml_content = {
            'pin_geometry': {'fuel_radius': 0.4, 'gap_radius': 0.41, 'clad_radius': 0.475},
            'lattice_description': [['UOX45']]
        }
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
        yaml.dump(yaml_content, f)
        f.close()
        return f.name

    @staticmethod
    def _create_minimal_scheme_yaml():
        yaml_content = {
            'calculation_steps': [
                {
                    'step_name': 'SSH',
                    'step_type': 'self_shielding',
                    'spatial_method': 'IC'
                }
            ]
        }
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
        yaml.dump(yaml_content, f)
        f.close()
        return f.name

    @staticmethod
    def _cleanup_yaml_files(*files):
        for f in files:
            try:
                os.unlink(f)
            except:
                pass


class TestDeprecationWarning:
    """Test deprecation warning for associate_temperatures_from_materials_yaml."""

    def test_deprecation_warning_emitted(self):
        """Verify deprecation warning is emitted"""
        mats_yaml = self._create_minimal_mats_yaml()
        geom_yaml = self._create_minimal_geom_yaml()

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = associate_temperatures_from_materials_yaml(mats_yaml, geom_yaml)

                # Check that at least one DeprecationWarning was raised
                assert len(w) >= 1
                assert issubclass(w[-1].category, DeprecationWarning)
                assert "deprecated" in str(w[-1].message).lower()
        finally:
            self._cleanup_yaml_files(mats_yaml, geom_yaml)

    @staticmethod
    def _create_minimal_mats_yaml():
        yaml_content = {
            'MIX_COMPOSITIONS': [
                {
                    'name': 'UOX',
                    'rod_id': 'UOX45',
                    'isotopic_composition': {'U235': 1e-4},
                    'temperature': 900.0  # Add temperature so function doesn't fail
                }
            ]
        }
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
        yaml.dump(yaml_content, f)
        f.close()
        return f.name

    @staticmethod
    def _create_minimal_geom_yaml():
        yaml_content = {
            'lattice_description': [['UOX45']]
        }
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
        yaml.dump(yaml_content, f)
        f.close()
        return f.name

    @staticmethod
    def _cleanup_yaml_files(*files):
        for f in files:
            try:
                os.unlink(f)
            except:
                pass


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
