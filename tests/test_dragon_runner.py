# Tests for the Dragon case execution runner.
# Covers: executable resolution, draglib resolution, scheme descriptor
#         building, run directory structure, manifest generation,
#         input staging, config archival, dry_run mode, latest symlink,
#         and keff/normal-end parsing.
#
# R.Guasch — 12/03/2026

import os
import textwrap

import pytest
import yaml

from starterDD.InterfaceToDD.dragon_runner import (
    DragonRunner,
    RunResult,
    _build_scheme_descriptor,
    _parse_keff_from_listing,
    _check_normal_end,
)
from starterDD.InterfaceToDD.case_generator import DragonCase
from starterDD.DDModel.DragonCalculationScheme import (
    CalculationStep,
    CalculationBranch,
    DragonCalculationScheme,
)
from conftest import (
    GE14_COMPOSITIONS_YAML,
    GE14_DOM_GEOMETRY_YAML,
    GE14_CALC_SCHEME_YAML,
    GE14_CALC_SCHEME_1L_YAML,
    GE14_TDT_DIR,
)


# =========================================================
# Scheme descriptor tests
# =========================================================

class TestSchemeDescriptor:
    """Tests for _build_scheme_descriptor."""

    def test_1L_MOC_scheme(self):
        """1-level MOC with RSE+IC self-shielding → RSE_IC_1L_MOC."""
        scheme = DragonCalculationScheme(name="test")
        scheme.add_step(CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS",
            self_shielding_method="RSE",
            spatial_method="IC",
        ))
        scheme.add_step(CalculationStep(
            name="FLUX", step_type="flux",
            spatial_method="MOC",
            tracking="TSPC",
            polar_angles_quadrature="GAUS",
            number_of_polar_angles=4,
        ))
        desc = _build_scheme_descriptor(scheme)
        assert desc == "RSE_IC_1L_MOC"

    def test_2L_IC_MOC_scheme(self):
        """2-level IC+MOC → RSE_IC_2L_IC_MOC."""
        scheme = DragonCalculationScheme(name="test")
        scheme.add_step(CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS",
            self_shielding_method="RSE",
            spatial_method="IC",
        ))
        scheme.add_step(CalculationStep(
            name="FLUX_L1", step_type="flux",
            spatial_method="IC", flux_level=1,
        ))
        scheme.add_step(CalculationStep(
            name="FLUX_L2", step_type="flux",
            spatial_method="MOC", flux_level=2,
            tracking="TSPC",
            polar_angles_quadrature="GAUS",
            number_of_polar_angles=4,
        ))
        desc = _build_scheme_descriptor(scheme)
        assert desc == "RSE_IC_2L_IC_MOC"

    def test_PT_CP_scheme(self):
        """PT+CP self-shielding, CP flux → PT_CP_1L_CP."""
        scheme = DragonCalculationScheme(name="test")
        scheme.add_step(CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS",
            self_shielding_method="PT",
            spatial_method="CP",
        ))
        scheme.add_step(CalculationStep(
            name="FLUX", step_type="flux",
            spatial_method="CP",
        ))
        desc = _build_scheme_descriptor(scheme)
        assert desc == "PT_CP_1L_CP"

    def test_with_statepoints(self):
        """Statepoint count appended when > 1."""
        scheme = DragonCalculationScheme(name="test")
        scheme.add_step(CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS",
            self_shielding_method="RSE",
            spatial_method="IC",
        ))
        scheme.add_step(CalculationStep(
            name="FLUX", step_type="flux",
            spatial_method="MOC",
            tracking="TSPC",
            polar_angles_quadrature="GAUS",
            number_of_polar_angles=4,
        ))
        scheme.branches = [
            CalculationBranch("CD", "coolant_density", [0.7, 0.5, 0.3]),
            CalculationBranch("FT", "fuel_temperature", [900.0]),
        ]
        desc = _build_scheme_descriptor(scheme)
        assert desc == "RSE_IC_1L_MOC_3sp"

    def test_single_statepoint_no_suffix(self):
        """Single statepoint (sp=1) → no suffix."""
        scheme = DragonCalculationScheme(name="test")
        scheme.add_step(CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS",
            self_shielding_method="RSE",
            spatial_method="IC",
        ))
        scheme.add_step(CalculationStep(
            name="FLUX", step_type="flux",
            spatial_method="MOC",
            tracking="TSPC",
            polar_angles_quadrature="GAUS",
            number_of_polar_angles=4,
        ))
        scheme.branches = [
            CalculationBranch("FT", "fuel_temperature", [900.0]),
        ]
        desc = _build_scheme_descriptor(scheme)
        assert desc == "RSE_IC_1L_MOC"

    def test_from_real_yaml(self):
        """Descriptor from the real CALC_SCHEME_1L.yaml file.
        SSH=RSE+IC, FLUX=MOC, 5 statepoints → RSE_IC_1L_MOC_5sp."""
        scheme = DragonCalculationScheme.from_yaml(
            GE14_CALC_SCHEME_1L_YAML
        )
        desc = _build_scheme_descriptor(scheme)
        assert desc == "RSE_IC_1L_MOC_5sp"


# =========================================================
# keff and normal-end parsing
# =========================================================

class TestListingParsing:
    """Tests for parsing Dragon listing output."""

    def test_parse_keff(self, tmp_path):
        """keff value is extracted from listing."""
        listing = tmp_path / "test.result"
        listing.write_text(
            "some output\n"
            "GE14_DOM keff= 1.08234\n"
            "more output\n"
        )
        assert _parse_keff_from_listing(str(listing)) == 1.08234

    def test_parse_keff_scientific(self, tmp_path):
        """keff in scientific notation."""
        listing = tmp_path / "test.result"
        listing.write_text("keff= 1.082e+00\n")
        assert abs(
            _parse_keff_from_listing(str(listing)) - 1.082
        ) < 1e-6

    def test_parse_keff_missing(self, tmp_path):
        """No keff in listing → None."""
        listing = tmp_path / "test.result"
        listing.write_text("no keff here\n")
        assert _parse_keff_from_listing(str(listing)) is None

    def test_parse_keff_nonexistent_file(self):
        """Non-existent file → None."""
        assert _parse_keff_from_listing("/nonexistent") is None

    def test_normal_end_success(self, tmp_path):
        """Listing with 'normal end' → True."""
        listing = tmp_path / "test.result"
        listing.write_text(
            "processing...\n" * 50 +
            " normal end of execution\n"
        )
        assert _check_normal_end(str(listing)) is True

    def test_normal_end_failure(self, tmp_path):
        """Listing without 'normal end' → False."""
        listing = tmp_path / "test.result"
        listing.write_text("processing...\nerror occurred\n")
        assert _check_normal_end(str(listing)) is False

    def test_normal_end_nonexistent(self):
        """Non-existent file → False."""
        assert _check_normal_end("/nonexistent") is False


# =========================================================
# Executable resolution
# =========================================================

class TestExecutableResolution:
    """Tests for DragonRunner._resolve_executable."""

    def test_explicit_path(self, tmp_path):
        """Explicit path to existing file is accepted."""
        exe = tmp_path / "dragon_exe"
        exe.write_text("#!/bin/sh\n")
        result = DragonRunner._resolve_executable(str(exe))
        assert result == str(exe.resolve())

    def test_env_var_dragon_exec(self, tmp_path, monkeypatch):
        """Resolves from $dragon_exec."""
        exe = tmp_path / "dragon_exe"
        exe.write_text("#!/bin/sh\n")
        monkeypatch.setenv("dragon_exec", str(exe))
        result = DragonRunner._resolve_executable(None)
        assert result == str(exe.resolve())

    def test_env_var_DRAGON_EXEC(self, tmp_path, monkeypatch):
        """Falls back to $DRAGON_EXEC."""
        exe = tmp_path / "dragon_exe"
        exe.write_text("#!/bin/sh\n")
        monkeypatch.delenv("dragon_exec", raising=False)
        monkeypatch.setenv("DRAGON_EXEC", str(exe))
        result = DragonRunner._resolve_executable(None)
        assert result == str(exe.resolve())

    def test_no_executable_raises(self, monkeypatch):
        """No valid path → FileNotFoundError."""
        monkeypatch.delenv("dragon_exec", raising=False)
        monkeypatch.delenv("DRAGON_EXEC", raising=False)
        with pytest.raises(FileNotFoundError, match="Cannot find"):
            DragonRunner._resolve_executable(None)

    def test_nonexistent_explicit_fallback_to_env(
        self, tmp_path, monkeypatch
    ):
        """Non-existent explicit path falls back to env."""
        exe = tmp_path / "real_dragon"
        exe.write_text("#!/bin/sh\n")
        monkeypatch.setenv("dragon_exec", str(exe))
        result = DragonRunner._resolve_executable(
            "/nonexistent/dragon"
        )
        assert result == str(exe.resolve())


# =========================================================
# Draglib resolution
# =========================================================

class TestDraglibResolution:
    """Tests for draglib path resolution.

    Since resolution is deferred to ``run()``, we test
    ``_resolve_draglib_paths`` directly.
    """

    def test_explicit_paths(self, tmp_path, monkeypatch):
        """Explicit draglib_paths dict works."""
        draglib = tmp_path / "draglibendfb8r1SHEM295"
        draglib.write_text("fake draglib")

        case = DragonCase(
            case_name="test",
            call_glow=False,
            draglib_name_to_alias={
                "draglibendfb8r1SHEM295": "endfb8r1_295",
            },
            config_yamls={
                "MATS": GE14_COMPOSITIONS_YAML,
                "GEOM": GE14_DOM_GEOMETRY_YAML,
                "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
            },
            output_path=str(tmp_path / "output"),
            tdt_path=GE14_TDT_DIR,
            tdt_base_name="GE14_DOM",
        )

        runner = DragonRunner(
            dragon_case=case,
            draglib_paths={
                "draglibendfb8r1SHEM295": str(draglib),
            },
        )
        resolved = runner._resolve_draglib_paths({
            "draglibendfb8r1SHEM295": str(draglib),
        })
        assert "draglibendfb8r1SHEM295" in resolved

    def test_draglibs_env_var(self, tmp_path, monkeypatch):
        """Resolves from $DRAGLIBS directory."""
        libs_dir = tmp_path / "libs"
        libs_dir.mkdir()
        draglib = libs_dir / "draglibendfb8r1SHEM295_v5p1"
        draglib.write_text("fake draglib")
        monkeypatch.delenv("DRAGLIBS_DIR", raising=False)
        monkeypatch.delenv("DRAGLIB_DIR", raising=False)
        monkeypatch.setenv("DRAGLIBS", str(libs_dir))

        case = DragonCase(
            case_name="test",
            call_glow=False,
            draglib_name_to_alias={
                "draglibendfb8r1SHEM295_v5p1": "endfb8r1_295",
            },
            config_yamls={
                "MATS": GE14_COMPOSITIONS_YAML,
                "GEOM": GE14_DOM_GEOMETRY_YAML,
                "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
            },
            output_path=str(tmp_path / "output"),
            tdt_path=GE14_TDT_DIR,
            tdt_base_name="GE14_DOM",
        )

        runner = DragonRunner(dragon_case=case)
        resolved = runner._resolve_draglib_paths(None)
        assert resolved[
            "draglibendfb8r1SHEM295_v5p1"
        ] == str(draglib.resolve())

    def test_missing_draglib_raises(self, tmp_path, monkeypatch):
        """Missing draglib → FileNotFoundError."""
        monkeypatch.delenv("DRAGLIBS_DIR", raising=False)
        monkeypatch.delenv("DRAGLIB_DIR", raising=False)
        monkeypatch.delenv("DRAGLIBS", raising=False)

        case = DragonCase(
            case_name="test",
            call_glow=False,
            draglib_name_to_alias={
                "draglibendfb8r1SHEM295_v5p1": "endfb8r1_295",
            },
            config_yamls={
                "MATS": GE14_COMPOSITIONS_YAML,
                "GEOM": GE14_DOM_GEOMETRY_YAML,
                "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
            },
            output_path=str(tmp_path / "output"),
            tdt_path=GE14_TDT_DIR,
            tdt_base_name="GE14_DOM",
        )

        runner = DragonRunner(dragon_case=case)
        with pytest.raises(FileNotFoundError, match="draglib"):
            runner._resolve_draglib_paths(None)

    def test_gz_draglib_resolved(self, tmp_path, monkeypatch):
        """A .gz compressed draglib is resolved."""
        libs_dir = tmp_path / "libs"
        libs_dir.mkdir()
        draglib_gz = libs_dir / "draglibendfb8r1SHEM295_v5p1.gz"
        draglib_gz.write_text("fake compressed")
        # Clear higher-priority env vars so DRAGLIBS is checked
        monkeypatch.delenv("DRAGLIBS_DIR", raising=False)
        monkeypatch.delenv("DRAGLIB_DIR", raising=False)
        monkeypatch.setenv("DRAGLIBS", str(libs_dir))

        case = DragonCase(
            case_name="test",
            call_glow=False,
            draglib_name_to_alias={
                "draglibendfb8r1SHEM295_v5p1": "endfb8r1_295",
            },
            config_yamls={
                "MATS": GE14_COMPOSITIONS_YAML,
                "GEOM": GE14_DOM_GEOMETRY_YAML,
                "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
            },
            output_path=str(tmp_path / "output"),
            tdt_path=GE14_TDT_DIR,
            tdt_base_name="GE14_DOM",
        )

        runner = DragonRunner(dragon_case=case)
        resolved = runner._resolve_draglib_paths(None)
        assert resolved[
            "draglibendfb8r1SHEM295_v5p1"
        ].endswith(".gz")


# =========================================================
# Run directory structure
# =========================================================

class TestRunDirectory:
    """Tests for run directory creation and structure."""

    def _make_runner(self, tmp_path, monkeypatch):
        """Create a DragonRunner with fake executable and draglib."""
        exe = tmp_path / "dragon_exe"
        exe.write_text("#!/bin/sh\n")
        draglib = tmp_path / "draglibendfb8r1SHEM295_v5p1"
        draglib.write_text("fake")
        monkeypatch.setenv("dragon_exec", str(exe))

        case = DragonCase(
            case_name="GE14_DOM",
            call_glow=False,
            draglib_name_to_alias={
                "draglibendfb8r1SHEM295_v5p1": "endfb8r1_295",
            },
            config_yamls={
                "MATS": GE14_COMPOSITIONS_YAML,
                "GEOM": GE14_DOM_GEOMETRY_YAML,
                "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
            },
            output_path=str(tmp_path / "output"),
            tdt_path=GE14_TDT_DIR,
            tdt_base_name="GE14_DOM",
        )
        # Generate procedures to populate scheme
        case.generate_cle2000_procedures()

        runner = DragonRunner(
            dragon_case=case,
            draglib_paths={
                "draglibendfb8r1SHEM295": str(draglib),
            },
            results_root=str(tmp_path / "results"),
        )
        runner._scheme = case.scheme
        runner._procedure_files = {
            "x2m": os.path.join(
                str(tmp_path / "output"), "GE14_DOM.x2m"
            ),
            "mix": os.path.join(
                str(tmp_path / "output"), "MIX_GE14_DOM.c2m"
            ),
            "trk": os.path.join(
                str(tmp_path / "output"), "TRK_GE14_DOM.c2m"
            ),
            "edir": os.path.join(
                str(tmp_path / "output"), "EDIR_GE14_DOM.c2m"
            ),
        }
        return runner

    def test_directory_structure(self, tmp_path, monkeypatch):
        """Run directory follows
        <root>/<case>/<scheme_desc>/<timestamp>/ structure."""
        runner = self._make_runner(tmp_path, monkeypatch)
        run_dir = runner._build_run_directory(
            "2026-03-12T14-32-07"
        )

        assert os.path.isdir(run_dir)
        assert "GE14_DOM" in run_dir
        assert "RSE_IC_1L_MOC_5sp" in run_dir
        assert "2026-03-12T14-32-07" in run_dir
        assert os.path.isdir(
            os.path.join(run_dir, "inputs")
        )
        assert os.path.isdir(
            os.path.join(run_dir, "procedures")
        )

    def test_config_yaml_archival(self, tmp_path, monkeypatch):
        """Config YAMLs are copied to inputs/ sub-directory."""
        runner = self._make_runner(tmp_path, monkeypatch)
        run_dir = runner._build_run_directory(
            "2026-03-12T14-32-07"
        )
        runner._archive_config_yamls(run_dir)

        inputs_dir = os.path.join(run_dir, "inputs")
        archived = os.listdir(inputs_dir)
        assert "material_compositions.yaml" in archived
        assert "GEOM_GE14_DOM.yaml" in archived
        assert "CALC_SCHEME_1L.yaml" in archived

    def test_latest_symlink(self, tmp_path, monkeypatch):
        """latest symlink points to most recent run directory."""
        runner = self._make_runner(tmp_path, monkeypatch)
        run_dir1 = runner._build_run_directory(
            "2026-03-12T14-00-00"
        )
        runner._update_latest_symlink(run_dir1)

        parent = os.path.dirname(run_dir1)
        latest = os.path.join(parent, "latest")
        assert os.path.islink(latest)
        assert os.readlink(latest) == "2026-03-12T14-00-00"

        # Second run updates the symlink
        run_dir2 = runner._build_run_directory(
            "2026-03-12T15-00-00"
        )
        runner._update_latest_symlink(run_dir2)
        assert os.readlink(latest) == "2026-03-12T15-00-00"


# =========================================================
# Manifest generation
# =========================================================

class TestManifest:
    """Tests for run manifest content."""

    def _make_runner(self, tmp_path, monkeypatch):
        """Create a DragonRunner with fake executable and draglib."""
        exe = tmp_path / "dragon_exe"
        exe.write_text("#!/bin/sh\n")
        draglib = tmp_path / "draglibendfb8r1SHEM295"
        draglib.write_text("fake")
        monkeypatch.setenv("dragon_exec", str(exe))

        case = DragonCase(
            case_name="GE14_DOM",
            call_glow=False,
            draglib_name_to_alias={
                "draglibendfb8r1SHEM295": "endfb8r1_295",
            },
            config_yamls={
                "MATS": GE14_COMPOSITIONS_YAML,
                "GEOM": GE14_DOM_GEOMETRY_YAML,
                "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
            },
            output_path=str(tmp_path / "output"),
            tdt_path=GE14_TDT_DIR,
            tdt_base_name="GE14_DOM",
        )
        case.generate_cle2000_procedures()

        runner = DragonRunner(
            dragon_case=case,
            draglib_paths={
                "draglibendfb8r1SHEM295": str(draglib),
            },
            results_root=str(tmp_path / "results"),
        )
        runner._scheme = case.scheme
        runner._procedure_files = {
            "x2m": os.path.join(
                str(tmp_path / "output"), "GE14_DOM.x2m"
            ),
        }
        return runner

    def test_manifest_content(self, tmp_path, monkeypatch):
        """Manifest contains all expected fields."""
        runner = self._make_runner(tmp_path, monkeypatch)
        run_dir = runner._build_run_directory(
            "2026-03-12T14-32-07"
        )
        manifest, manifest_path = runner._write_run_manifest(
            run_dir, "2026-03-12T14-32-07"
        )

        assert os.path.isfile(manifest_path)
        assert manifest["case_name"] == "GE14_DOM"
        assert manifest["assembly_id"] == "GE14_DOM"
        assert manifest["run_timestamp"] == "2026-03-12T14-32-07"
        assert manifest["num_threads"] == 1
        assert manifest["execution"] is None

        # Scheme info
        scheme_info = manifest["calculation_scheme"]
        assert scheme_info["descriptor"] == "RSE_IC_1L_MOC_5sp"
        assert len(scheme_info["steps"]) == 2
        assert scheme_info["steps"][0]["ssh_method"] == "RSE"
        assert scheme_info["steps"][0]["method"] == "IC"
        assert scheme_info["steps"][1]["method"] == "MOC"
        assert scheme_info["total_statepoints"] == 5
        assert "coolant_density" in scheme_info["branches"]

        # Config yaml references
        assert "MATS" in manifest["config_yamls"]
        assert "GEOM" in manifest["config_yamls"]
        assert "CALC_SCHEME" in manifest["config_yamls"]

    def test_manifest_yaml_roundtrip(self, tmp_path, monkeypatch):
        """Manifest file can be loaded back as YAML."""
        runner = self._make_runner(tmp_path, monkeypatch)
        run_dir = runner._build_run_directory("2026-03-12T14-32-07")
        _, manifest_path = runner._write_run_manifest(
            run_dir, "2026-03-12T14-32-07"
        )

        with open(manifest_path) as f:
            loaded = yaml.safe_load(f)

        assert loaded["case_name"] == "GE14_DOM"
        assert loaded["calculation_scheme"]["descriptor"] == (
            "RSE_IC_1L_MOC_5sp"
        )

    def test_manifest_updated_with_results(
        self, tmp_path, monkeypatch
    ):
        """After execution, manifest is updated with results."""
        runner = self._make_runner(tmp_path, monkeypatch)
        run_dir = runner._build_run_directory("2026-03-12T14-32-07")
        manifest, manifest_path = runner._write_run_manifest(
            run_dir, "2026-03-12T14-32-07"
        )

        result = RunResult(
            success=True,
            run_directory=run_dir,
            manifest_path=manifest_path,
            wall_time_seconds=127.4,
            return_code=0,
            listing_path="test.result",
            keff=1.08234,
        )
        runner._update_manifest_with_results(
            manifest, manifest_path, result
        )

        with open(manifest_path) as f:
            loaded = yaml.safe_load(f)

        assert loaded["execution"]["success"] is True
        assert loaded["execution"]["wall_time_seconds"] == 127.4
        assert loaded["execution"]["keff"] == 1.08234
        assert loaded["execution"]["return_code"] == 0


# =========================================================
# Input staging
# =========================================================

class TestInputStaging:
    """Tests for input file staging into execution directory."""

    def _make_runner_with_procs(self, tmp_path, monkeypatch):
        """Create runner with real generated procedures."""
        draglib = tmp_path / "draglibendfb8r1SHEM295"
        draglib.write_text("fake draglib content")

        output_dir = tmp_path / "output"
        case = DragonCase(
            case_name="GE14_DOM",
            call_glow=False,
            draglib_name_to_alias={
                "draglibendfb8r1SHEM295": "endfb8r1_295",
            },
            config_yamls={
                "MATS": GE14_COMPOSITIONS_YAML,
                "GEOM": GE14_DOM_GEOMETRY_YAML,
                "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
            },
            output_path=str(output_dir),
            tdt_path=GE14_TDT_DIR,
            tdt_base_name="GE14_DOM",
        )
        procs = case.generate_cle2000_procedures()

        runner = DragonRunner(
            dragon_case=case,
            draglib_paths={
                "draglibendfb8r1SHEM295": str(draglib),
            },
            results_root=str(tmp_path / "results"),
        )
        runner._scheme = case.scheme
        runner._procedure_files = procs
        # Resolve draglib paths for staging tests
        runner.draglib_paths = runner._resolve_draglib_paths({
            "draglibendfb8r1SHEM295": str(draglib),
        })
        return runner, case

    def test_procedures_copied_to_staging(
        self, tmp_path, monkeypatch
    ):
        """CLE2000 files are copied into the staging directory."""
        runner, _ = self._make_runner_with_procs(
            tmp_path, monkeypatch
        )
        run_dir = runner._build_run_directory("2026-03-12T14-32-07")
        staging_dir = str(tmp_path / "staging")
        os.makedirs(staging_dir)

        runner._stage_inputs(run_dir, staging_dir)

        staged_files = os.listdir(staging_dir)
        assert "GE14_DOM.x2m" in staged_files
        assert "MIX.c2m" in staged_files
        assert "TRK.c2m" in staged_files
        assert "EDIR.c2m" in staged_files

    def test_procedures_archived_in_run_dir(
        self, tmp_path, monkeypatch
    ):
        """CLE2000 files are also archived in run_dir/procedures/."""
        runner, _ = self._make_runner_with_procs(
            tmp_path, monkeypatch
        )
        run_dir = runner._build_run_directory("2026-03-12T14-32-07")
        staging_dir = str(tmp_path / "staging")
        os.makedirs(staging_dir)

        runner._stage_inputs(run_dir, staging_dir)

        proc_dir = os.path.join(run_dir, "procedures")
        archived = os.listdir(proc_dir)
        assert "GE14_DOM.x2m" in archived

    def test_tdt_files_symlinked(self, tmp_path, monkeypatch):
        """TDT .dat files are symlinked into staging directory."""
        runner, _ = self._make_runner_with_procs(
            tmp_path, monkeypatch
        )
        run_dir = runner._build_run_directory("2026-03-12T14-32-07")
        staging_dir = str(tmp_path / "staging")
        os.makedirs(staging_dir)

        runner._stage_inputs(run_dir, staging_dir)

        staged_files = os.listdir(staging_dir)
        dat_files = [f for f in staged_files if f.endswith(".dat")]
        assert len(dat_files) > 0
        # Should be symlinks
        for f in dat_files:
            assert os.path.islink(
                os.path.join(staging_dir, f)
            )

    def test_draglib_symlinked(self, tmp_path, monkeypatch):
        """Draglib file is symlinked into staging directory."""
        runner, _ = self._make_runner_with_procs(
            tmp_path, monkeypatch
        )
        run_dir = runner._build_run_directory("2026-03-12T14-32-07")
        staging_dir = str(tmp_path / "staging")
        os.makedirs(staging_dir)

        runner._stage_inputs(run_dir, staging_dir)

        staged_files = os.listdir(staging_dir)
        assert "endfb8r1_295" in staged_files
        assert os.path.islink(
            os.path.join(staging_dir, "endfb8r1_295")
        )


# =========================================================
# Dry run mode
# =========================================================

class TestDryRun:
    """Tests for dry_run=True mode."""

    def test_dry_run_creates_directory_no_execution(
        self, tmp_path, monkeypatch
    ):
        """dry_run stages everything but does not execute Dragon.
        
        No real executable or draglib files are needed for a
        dry run — this is the whole point of deferred resolution.
        """
        monkeypatch.delenv("dragon_exec", raising=False)
        monkeypatch.delenv("DRAGON_EXEC", raising=False)
        monkeypatch.delenv("DRAGLIBS", raising=False)

        output_dir = tmp_path / "output"
        case = DragonCase(
            case_name="GE14_DOM",
            call_glow=False,
            draglib_name_to_alias={
                "draglibendfb8r1SHEM295": "endfb8r1_295",
            },
            config_yamls={
                "MATS": GE14_COMPOSITIONS_YAML,
                "GEOM": GE14_DOM_GEOMETRY_YAML,
                "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
            },
            output_path=str(output_dir),
            tdt_path=GE14_TDT_DIR,
            tdt_base_name="GE14_DOM",
        )

        runner = DragonRunner(
            dragon_case=case,
            results_root=str(tmp_path / "results"),
        )

        result = runner.run(dry_run=True)

        assert result.success is True
        assert result.wall_time_seconds == 0.0
        assert result.return_code == 0
        assert os.path.isdir(result.run_directory)

        # Manifest should exist with dry_run flag
        with open(result.manifest_path) as f:
            manifest = yaml.safe_load(f)
        assert manifest["execution"]["dry_run"] is True

        # Config YAMLs should be archived
        inputs_dir = os.path.join(
            result.run_directory, "inputs"
        )
        assert os.path.isdir(inputs_dir)
        assert len(os.listdir(inputs_dir)) == 3

        # Procedures should be archived
        proc_dir = os.path.join(
            result.run_directory, "procedures"
        )
        assert os.path.isdir(proc_dir)
        assert "GE14_DOM.x2m" in os.listdir(proc_dir)

        # Staging contents file should exist
        staging_listing = os.path.join(
            result.run_directory, "_staging_contents.txt"
        )
        assert os.path.isfile(staging_listing)


# =========================================================
# DragonCase.run() convenience method
# =========================================================

class TestDragonCaseRun:
    """Tests for the DragonCase.run() convenience method."""

    def test_run_dry_run(self, tmp_path, monkeypatch):
        """DragonCase.run(dry_run=True) generates and stages
        without needing a real executable or draglibs."""
        monkeypatch.delenv("dragon_exec", raising=False)
        monkeypatch.delenv("DRAGON_EXEC", raising=False)
        monkeypatch.delenv("DRAGLIBS", raising=False)

        case = DragonCase(
            case_name="GE14_DOM",
            call_glow=False,
            draglib_name_to_alias={
                "draglibendfb8r1SHEM295": "endfb8r1_295",
            },
            config_yamls={
                "MATS": GE14_COMPOSITIONS_YAML,
                "GEOM": GE14_DOM_GEOMETRY_YAML,
                "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
            },
            output_path=str(tmp_path / "output"),
            tdt_path=GE14_TDT_DIR,
            tdt_base_name="GE14_DOM",
        )

        result = case.run(
            results_root=str(tmp_path / "results"),
            dry_run=True,
        )

        assert result.success is True
        assert os.path.isdir(result.run_directory)
        assert "GE14_DOM" in result.run_directory
        assert "RSE_IC_1L_MOC_5sp" in result.run_directory


# =========================================================
# Output collection
# =========================================================

class TestOutputCollection:
    """Tests for _collect_outputs."""

    def test_result_file_moved(self, tmp_path):
        """The .result listing is moved to run_dir."""
        staging = tmp_path / "staging"
        staging.mkdir()
        run_dir = tmp_path / "run"
        run_dir.mkdir()

        listing = staging / "GE14_DOM.result"
        listing.write_text("test output")

        # Create a minimal runner-like object with case_name
        class FakeCase:
            case_name = "GE14_DOM"

        class FakeRunner(DragonRunner):
            def __init__(self):
                self.case = FakeCase()

        runner = FakeRunner.__new__(FakeRunner)
        runner.case = FakeCase()
        runner._collect_outputs(str(staging), str(run_dir))

        assert os.path.isfile(
            str(run_dir / "GE14_DOM.result")
        )
        assert not os.path.isfile(str(listing))

    def test_compo_files_moved(self, tmp_path):
        """_CPO_* files are moved to run_dir."""
        staging = tmp_path / "staging"
        staging.mkdir()
        run_dir = tmp_path / "run"
        run_dir.mkdir()

        compo = staging / "_CPO_GE14_DOM"
        compo.write_text("compo data")

        class FakeCase:
            case_name = "GE14_DOM"

        runner = DragonRunner.__new__(DragonRunner)
        runner.case = FakeCase()
        runner._collect_outputs(str(staging), str(run_dir))

        assert os.path.isfile(
            str(run_dir / "_CPO_GE14_DOM")
        )

    def test_ps_files_moved(self, tmp_path):
        """PostScript files are moved to run_dir."""
        staging = tmp_path / "staging"
        staging.mkdir()
        run_dir = tmp_path / "run"
        run_dir.mkdir()

        ps = staging / "geometry.ps"
        ps.write_text("PS data")

        class FakeCase:
            case_name = "GE14_DOM"

        runner = DragonRunner.__new__(DragonRunner)
        runner.case = FakeCase()
        runner._collect_outputs(str(staging), str(run_dir))

        assert os.path.isfile(str(run_dir / "geometry.ps"))
