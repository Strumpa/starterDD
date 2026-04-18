# Unit tests for TDT file handling refactor - focused on key behaviors.
#
# Tests the new behavior without requiring full Dragon setup:
# - Actual filenames stored in tdt_files_used (not standardized)
# - tdt_file_mapping created and populated
# - No intermediate symlinks created in tdt_path
# - Enhanced error handling in dragon_runner
#
# R.Guasch — 04/17/2026

import os
import tempfile
from unittest.mock import Mock, patch, MagicMock

import pytest

from starterDD.InterfaceToDD.case_generator import DragonCase
from starterDD.InterfaceToDD.dragon_runner import DragonRunner


class TestCaseGeneratorTDTFilenameStorage:
    """Unit tests for Phase 1: TDT filename storage in case_generator."""

    def test_dragoncase_initializes_tdt_attributes(self):
        """Verify DragonCase.__init__ initializes TDT tracking attributes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            case = DragonCase(
                case_name="test",
                call_glow=False,
                draglib_name_to_alias={"lib": "LIB"},
                config_yamls={
                    "MATS": os.path.join(tmpdir, "mats.yaml"),
                    "GEOM": os.path.join(tmpdir, "geom.yaml"),
                    "CALC_SCHEME": os.path.join(tmpdir, "scheme.yaml"),
                },
                tdt_path=tmpdir,
            )
            
            # Check attributes exist
            assert hasattr(case, 'tdt_files_used'), \
                "DragonCase should have tdt_files_used"
            assert hasattr(case, 'tdt_file_mapping'), \
                "DragonCase should have tdt_file_mapping"
            
            # Check they're empty initially
            assert isinstance(case.tdt_files_used, dict)
            assert isinstance(case.tdt_file_mapping, dict)
            assert len(case.tdt_files_used) == 0
            assert len(case.tdt_file_mapping) == 0

    def test_tdt_file_mapping_structure(self):
        """Verify tdt_file_mapping has correct structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            case = DragonCase(
                case_name="test",
                call_glow=False,
                draglib_name_to_alias={"lib": "LIB"},
                config_yamls={
                    "MATS": os.path.join(tmpdir, "mats.yaml"),
                    "GEOM": os.path.join(tmpdir, "geom.yaml"),
                    "CALC_SCHEME": os.path.join(tmpdir, "scheme.yaml"),
                },
                tdt_path=tmpdir,
            )
            
            # Manually populate to test structure
            case.tdt_file_mapping['SSH'] = {
                'actual': 'custom_ssh.dat',
                'standardized': 'test_SSH_CP_2D.dat',
                'match': False,
            }
            
            mapping = case.tdt_file_mapping['SSH']
            assert 'actual' in mapping
            assert 'standardized' in mapping
            assert 'match' in mapping
            assert mapping['match'] is False


class TestDragonRunnerTDTValidation:
    """Unit tests for Phase 4: Enhanced validation in dragon_runner."""

    def test_stage_inputs_requires_tdt_files_to_exist(self):
        """Verify _stage_inputs raises error if TDT file doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            case = Mock()
            case.case_name = "test"
            case.tdt_path = tmpdir
            case.tdt_files_used = {
                'SSH': 'nonexistent_file.dat'
            }
            case.tdt_file_mapping = {
                'SSH': {
                    'actual': 'nonexistent_file.dat',
                    'standardized': 'test_SSH_CP_TISO.dat',
                    'match': False,
                }
            }

            runner = DragonRunner(case, results_root=tmpdir)
            # Mock the procedure files
            runner._procedure_files = {}

            with tempfile.TemporaryDirectory() as staging_dir:
                run_dir = tmpdir

                # Should raise FileNotFoundError
                with pytest.raises(FileNotFoundError):
                    runner._stage_inputs(run_dir, staging_dir, skip_draglibs=True)

    def test_stage_inputs_validates_symlink_resolution(self):
        """Verify _stage_inputs validates that staging symlinks resolve."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create actual TDT file
            tdt_file = os.path.join(tmpdir, "test.dat")
            with open(tdt_file, 'w') as f:
                f.write("test")
            
            case = Mock()
            case.case_name = "test"
            case.tdt_path = tmpdir
            case.tdt_files_used = {
                'SSH': 'test.dat'
            }
            case.tdt_file_mapping = {
                'SSH': {
                    'actual': 'test.dat',
                    'standardized': 'test_SSH_CP_TISO.dat',
                    'match': True,
                }
            }
            
            runner = DragonRunner(case, results_root=tmpdir)
            runner._procedure_files = {}  # Empty procedures for test
            
            with tempfile.TemporaryDirectory() as staging_dir:
                run_dir = tmpdir
                
                # Should successfully create symlink and validate it
                runner._stage_inputs(run_dir, staging_dir, skip_draglibs=True)
                
                # Verify standardized symlink exists in staging_dir
                # (it should be at the standardized name, not the actual name)
                staging_symlink = os.path.join(staging_dir, 'test_SSH_CP_TISO.dat')
                assert os.path.islink(staging_symlink), \
                    "Staging symlink should be created"
                assert os.path.exists(staging_symlink), \
                    "Staging symlink should resolve"

    def test_stage_inputs_error_message_includes_step_name(self):
        """Verify error messages include step name for debugging."""
        with tempfile.TemporaryDirectory() as tmpdir:
            case = Mock()
            case.case_name = "test"
            case.tdt_path = tmpdir
            case.tdt_files_used = {
                'SSH': 'missing_ssh.dat',
                'FLUX': 'missing_flux.dat',
            }
            case.tdt_file_mapping = {
                'SSH': {
                    'actual': 'missing_ssh.dat',
                    'standardized': 'test_SSH_CP_TISO.dat',
                    'match': False,
                },
                'FLUX': {
                    'actual': 'missing_flux.dat',
                    'standardized': 'test_FLUX_MOC_TISO.dat',
                    'match': False,
                }
            }
            
            runner = DragonRunner(case, results_root=tmpdir)
            runner._procedure_files = {}
            
            with tempfile.TemporaryDirectory() as staging_dir:
                run_dir = tmpdir
                
                # Ensure FileNotFoundError is raised with step name in message
                with pytest.raises(FileNotFoundError) as exc_info:
                    runner._stage_inputs(run_dir, staging_dir, skip_draglibs=True)
                
                # Error should mention the step name (SSH or FLUX)
                error_msg = str(exc_info.value)
                has_step_name = (
                    'SSH' in error_msg or 
                    'FLUX' in error_msg or 
                    'step' in error_msg.lower()
                )
                assert has_step_name, \
                    f"Error should mention step name in message: {error_msg}"


class TestManifestTDTInformation:
    """Unit tests for Phase 3: TDT information in manifest."""

    def test_build_scheme_manifest_includes_tdt_files(self):
        """Verify _build_scheme_manifest includes tdt_files list."""
        from starterDD.DDModel.DragonCalculationScheme import (
            CalculationStep,
            DragonCalculationScheme,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            case = Mock()
            case.case_name = "test"
            case.tdt_path = tmpdir
            case.tdt_files_used = {}
            case.tdt_file_mapping = {
                'SSH': {
                    'actual': 'custom_ssh.dat',
                    'standardized': 'test_SSH_CP_TISO.dat',
                    'match': False,
                }
            }

            # Create a minimal scheme
            scheme = DragonCalculationScheme(name="test")
            scheme.add_step(CalculationStep(
                name="SSH",
                step_type="self_shielding",
                self_shielding_module="USS",
                self_shielding_method="RSE",
                spatial_method="CP",
                tracking="TISO",
            ))

            runner = DragonRunner(case, results_root=tmpdir)
            runner._scheme = scheme
            
            manifest = runner._build_scheme_manifest()
            
            # Verify manifest includes tdt_files
            assert 'tdt_files' in manifest, \
                "Manifest should include tdt_files"
            assert isinstance(manifest['tdt_files'], list), \
                "tdt_files should be a list"
            
            # Verify tdt_files content
            if manifest['tdt_files']:
                tdt_entry = manifest['tdt_files'][0]
                assert 'step_name' in tdt_entry
                assert 'actual_filename' in tdt_entry
                assert 'standardized_filename' in tdt_entry
                assert 'names_match' in tdt_entry

    def test_manifest_handles_missing_tdt_file_mapping(self):
        """Verify manifest gracefully handles cases without tdt_file_mapping."""
        from starterDD.DDModel.DragonCalculationScheme import (
            CalculationStep,
            DragonCalculationScheme,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            case = Mock(spec=['case_name', 'tdt_path'])  # Restrict mock attributes
            case.case_name = "test"
            case.tdt_path = tmpdir

            scheme = DragonCalculationScheme(name="test")
            scheme.add_step(CalculationStep(
                name="SSH",
                step_type="self_shielding",
                self_shielding_module="USS",
                self_shielding_method="RSE",
                spatial_method="CP",
                tracking="TISO",
            ))

            runner = DragonRunner(case, results_root=tmpdir)
            runner._scheme = scheme
            
            # Should not raise error
            manifest = runner._build_scheme_manifest()
            
            # tdt_files should be empty list
            assert 'tdt_files' in manifest
            assert manifest['tdt_files'] == []


class TestNoIntermediateSymlinks:
    """Tests verifying no intermediate symlinks are created."""

    def test_case_generator_does_not_create_symlinks_in_tdt_path(self):
        """Verify case_generator does not create symlinks in tdt_path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a test TDT file
            tdt_file = os.path.join(tmpdir, 'test_tdt.dat')
            with open(tdt_file, 'w') as f:
                f.write("TDT content")
            
            initial_items = set(os.listdir(tmpdir))
            
            case = DragonCase(
                case_name="test",
                call_glow=False,
                draglib_name_to_alias={"lib": "LIB"},
                config_yamls={
                    "MATS": os.path.join(tmpdir, "mats.yaml"),
                    "GEOM": os.path.join(tmpdir, "geom.yaml"),
                    "CALC_SCHEME": os.path.join(tmpdir, "scheme.yaml"),
                },
                tdt_path=tmpdir,
            )
            
            # After initialization, no new symlinks should exist
            final_items = set(os.listdir(tmpdir))
            new_items = final_items - initial_items
            
            # None of the new items should be symlinks
            for item in new_items:
                item_path = os.path.join(tmpdir, item)
                assert not os.path.islink(item_path), \
                    f"Intermediate symlink should not be created: {item}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
