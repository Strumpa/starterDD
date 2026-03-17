## Tests for DragonCalculationScheme, CalculationStep, SectorConfig,
## BoxDiscretizationConfig and their integration with the assembly model.

import pytest
import tempfile
import os
import yaml

from starterDD.DDModel.DragonCalculationScheme import (
    DragonCalculationScheme,
    CalculationStep,
    SectorConfig,
    BoxDiscretizationConfig,
    CrossModeratorDiscretizationConfig,
    ControlCrossSubmeshConfig,
    VALID_STEP_TYPES,
    VALID_SELF_SHIELDING_MODULES,
    VALID_SELF_SHIELDING_METHODS,
    VALID_SPATIAL_METHODS,
    VALID_TRACKING_OPTIONS,
    VALID_RADIAL_SCHEMES,
)
from starterDD.DDModel.DragonModel import CartesianAssemblyModel, FuelPinModel
from starterDD.DDModel.helpers import associate_material_to_rod_ID
from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml

from conftest import (
    GE14_COMPOSITIONS_YAML,
    GE14_DOM_GEOMETRY_YAML,
    GE14_SIMPLE_GEOMETRY_YAML,
    GE14_CALC_SCHEME_YAML,
)


# =====================================================================
#  Helpers
# =====================================================================

# Paths anchored to the tests/ directory via conftest.py
_COMPOSITIONS_YAML = GE14_COMPOSITIONS_YAML
_GEOMETRY_DOM_YAML = GE14_DOM_GEOMETRY_YAML
_GEOMETRY_SIMPLE_YAML = GE14_SIMPLE_GEOMETRY_YAML
_CALC_SCHEME_YAML = GE14_CALC_SCHEME_YAML


def _build_simple_assembly():
    """
    Build and return a simplified GE14 assembly model with by_pin numbering
    applied (1 zone per pin, no Santamarina subdivision yet).
    """
    ROD_to_material = associate_material_to_rod_ID(_COMPOSITIONS_YAML,
                                                   _GEOMETRY_SIMPLE_YAML)
    assembly = CartesianAssemblyModel(
        name="simple_test",
        tdt_file="dummy.tdt",
        geometry_description_yaml=_GEOMETRY_SIMPLE_YAML,
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.set_uniform_temperatures(
        fuel_temperature=900.0, gap_temperature=600.0,
        coolant_temperature=600.0, moderator_temperature=600.0,
        structural_temperature=600.0,
    )
    assembly.analyze_lattice_description(build_pins=True)
    return assembly


# =====================================================================
#  SectorConfig
# =====================================================================

class TestSectorConfig:
    def test_defaults(self):
        sc = SectorConfig()
        assert sc.sectors == []
        assert sc.angles == []
        assert sc.windmill is False

    def test_custom(self):
        sc = SectorConfig(sectors=[1, 1, 8], angles=[0, 0, 22.5], windmill=True)
        assert sc.sectors == [1, 1, 8]
        assert sc.angles == [0, 0, 22.5]
        assert sc.windmill is True

    def test_repr(self):
        sc = SectorConfig(sectors=[4], angles=[0])
        r = repr(sc)
        assert "SectorConfig" in r
        assert "[4]" in r

    def test_splits_default_is_none(self):
        """SectorConfig.splits defaults to None when not provided."""
        sc = SectorConfig()
        assert sc.splits is None

    def test_splits_stored_as_tuple(self):
        """SectorConfig.splits should be stored as a tuple."""
        sc = SectorConfig(splits=[3, 3])
        assert sc.splits == (3, 3)
        assert isinstance(sc.splits, tuple)

    def test_splits_repr(self):
        """repr should show splits when set."""
        sc = SectorConfig(splits=[4, 5])
        r = repr(sc)
        assert "splits=(4, 5)" in r
        # Should NOT show sectors/angles in the splits repr
        assert "sectors=" not in r

    def test_splits_with_sectors_both_stored(self):
        """Both splits and sectors can be set simultaneously (warning at build time)."""
        sc = SectorConfig(sectors=[1, 1, 8], angles=[0, 0, 22.5], splits=[3, 3])
        assert sc.splits == (3, 3)
        assert sc.sectors == [1, 1, 8]

    # ----- additional_radial_splits_in_moderator -----

    def test_additional_radial_splits_default(self):
        """Default additional_radial_splits_in_moderator is 1 (no extra splits)."""
        sc = SectorConfig()
        assert sc.additional_radial_splits_in_moderator == None

    def test_additional_radial_splits_int(self):
        """Int value is stored as-is."""
        sc = SectorConfig(additional_radial_splits_in_moderator=3)
        assert sc.additional_radial_splits_in_moderator == 3

    def test_additional_radial_splits_list(self):
        """List of radii is stored as-is."""
        sc = SectorConfig(additional_radial_splits_in_moderator=[0.3, 0.6])
        assert sc.additional_radial_splits_in_moderator == [0.3, 0.6]

    # ----- resolve_water_rod_radii -----

    def test_resolve_radii_default_no_splits(self):
        """N=1 → empty list (no extra radii)."""
        sc = SectorConfig(additional_radial_splits_in_moderator=1)
        assert sc.resolve_water_rod_radii(1.0) == []

    def test_resolve_radii_none_no_splits(self):
        """None → empty list."""
        sc = SectorConfig(additional_radial_splits_in_moderator=None)
        assert sc.resolve_water_rod_radii(1.0) == []

    def test_resolve_radii_int_2(self):
        """N=2 → 1 intermediate radius at inner_radius/2."""
        sc = SectorConfig(additional_radial_splits_in_moderator=2)
        radii = sc.resolve_water_rod_radii(1.0)
        assert len(radii) == 1
        assert radii[0] == pytest.approx(0.5)

    def test_resolve_radii_int_3(self):
        """N=3 → 2 intermediate radii at 1/3 and 2/3 of inner_radius."""
        sc = SectorConfig(additional_radial_splits_in_moderator=3)
        radii = sc.resolve_water_rod_radii(1.2)
        assert len(radii) == 2
        assert radii[0] == pytest.approx(0.4)
        assert radii[1] == pytest.approx(0.8)

    def test_resolve_radii_user_defined_list(self):
        """Explicit radii list is returned sorted."""
        sc = SectorConfig(additional_radial_splits_in_moderator=[0.6, 0.3])
        radii = sc.resolve_water_rod_radii(1.0)
        assert radii == [0.3, 0.6]

    def test_resolve_radii_user_defined_invalid_too_large(self):
        """Radius >= inner_radius should raise ValueError."""
        sc = SectorConfig(additional_radial_splits_in_moderator=[0.5, 1.0])
        with pytest.raises(ValueError, match="out of range"):
            sc.resolve_water_rod_radii(1.0)

    def test_resolve_radii_user_defined_invalid_negative(self):
        """Radius <= 0 should raise ValueError."""
        sc = SectorConfig(additional_radial_splits_in_moderator=[-0.1, 0.5])
        with pytest.raises(ValueError, match="out of range"):
            sc.resolve_water_rod_radii(1.0)

    def test_resolve_radii_invalid_type(self):
        """Non-int, non-list type should raise TypeError."""
        sc = SectorConfig(additional_radial_splits_in_moderator="two")
        with pytest.raises(TypeError, match="must be int or list"):
            sc.resolve_water_rod_radii(1.0)

    # ----- expanded_sectors_and_angles -----

    def test_expanded_no_extra(self):
        """No extra radii → sectors and angles unchanged."""
        sc = SectorConfig(
            sectors=[16, 16, 16], angles=[0, 0, 0],
            additional_radial_splits_in_moderator=1,
        )
        s, a = sc.expanded_sectors_and_angles(1.0)
        assert s == [16, 16, 16]
        assert a == [0, 0, 0]

    def test_expanded_with_int_2(self):
        """N=2 → 1 extra ring; sectors[0]/angles[0] prepended once."""
        sc = SectorConfig(
            sectors=[16, 16, 16], angles=[0, 0, 0],
            additional_radial_splits_in_moderator=2,
        )
        s, a = sc.expanded_sectors_and_angles(1.0)
        assert s == [16, 16, 16, 16]
        assert a == [0, 0, 0, 0]

    def test_expanded_with_int_3(self):
        """N=3 → 2 extra rings; sectors[0]/angles[0] prepended twice."""
        sc = SectorConfig(
            sectors=[8, 4, 12], angles=[10.0, 20.0, 30.0],
            additional_radial_splits_in_moderator=3,
        )
        s, a = sc.expanded_sectors_and_angles(1.0)
        assert s == [8, 8, 8, 4, 12]
        assert a == [10.0, 10.0, 10.0, 20.0, 30.0]

    def test_expanded_with_user_defined_list(self):
        """Explicit 2-radii list → 2 extra entries prepended."""
        sc = SectorConfig(
            sectors=[16, 16, 16], angles=[0, 0, 0],
            additional_radial_splits_in_moderator=[0.3, 0.6],
        )
        s, a = sc.expanded_sectors_and_angles(1.0)
        assert s == [16, 16, 16, 16, 16]
        assert a == [0, 0, 0, 0, 0]

    def test_expanded_empty_sectors(self):
        """Empty sectors → returned unchanged regardless of splits."""
        sc = SectorConfig(
            sectors=[], angles=[],
            additional_radial_splits_in_moderator=3,
        )
        s, a = sc.expanded_sectors_and_angles(1.0)
        assert s == []
        assert a == []

    def test_repr_with_additional_splits(self):
        """repr should mention additional_radial_splits when != 1."""
        sc = SectorConfig(
            sectors=[16, 16, 16], angles=[0, 0, 0],
            additional_radial_splits_in_moderator=3,
        )
        r = repr(sc)
        assert "additional_radial_splits=3" in r

    def test_repr_without_additional_splits(self):
        """repr should NOT mention additional_radial_splits when == 1."""
        sc = SectorConfig(
            sectors=[16, 16, 16], angles=[0, 0, 0],
            additional_radial_splits_in_moderator=1,
        )
        r = repr(sc)
        assert "additional_radial_splits" not in r


# =====================================================================
#  BoxDiscretizationConfig
# =====================================================================

class TestBoxDiscretizationConfig:
    def test_defaults(self):
        bdc = BoxDiscretizationConfig()
        assert bdc.enabled is False
        assert bdc.corner_splits == (4, 4)
        assert bdc.gap_splits is None

    def test_custom(self):
        bdc = BoxDiscretizationConfig(
            enabled=True,
            corner_splits=[3, 3],
            gap_splits=[10, 1],
        )
        assert bdc.enabled is True
        assert bdc.corner_splits == (3, 3)
        assert bdc.gap_splits == (10, 1)

    def test_resolve_splits_defaults(self):
        bdc = BoxDiscretizationConfig(enabled=True)
        corner, side_h, side_v = bdc.resolve_splits(n_cols=10, n_rows=10)
        assert corner == (4, 4)
        assert side_h == (10, 1)
        assert side_v == (1, 10)

    def test_resolve_splits_explicit(self):
        bdc = BoxDiscretizationConfig(
            enabled=True,
            corner_splits=[2, 2],
            gap_splits=[5, 2],
        )
        corner, side_h, side_v = bdc.resolve_splits(n_cols=10, n_rows=10)
        assert corner == (2, 2)
        assert side_h == (5, 2)
        assert side_v == (2, 5)

    def test_deprecated_side_x_splits(self):
        """side_x_splits should be accepted as deprecated alias."""
        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bdc = BoxDiscretizationConfig(
                enabled=True,
                corner_splits=[3, 3],
                side_x_splits=[10, 1],
            )
        assert bdc.gap_splits == (10, 1)
        assert any(issubclass(warning.category, DeprecationWarning) for warning in w)


# =====================================================================
#  ControlCrossSubmeshConfig
# =====================================================================

class TestControlCrossSubmeshConfig:
    def test_defaults(self):
        cfg = ControlCrossSubmeshConfig()
        assert cfg.enabled is True
        assert cfg.control_cross_corner_splits is None
        assert cfg.central_structure_splits is None
        assert cfg.extend_splits_at_tube_boundaries is True
        assert cfg.split_tubes_in_half is False

    def test_custom(self):
        cfg = ControlCrossSubmeshConfig(
            enabled=True,
            control_cross_corner_splits=[2, 3],
            central_structure_splits=[4, 2],
            extend_splits_at_tube_boundaries=False,
            split_tubes_in_half=True,
        )
        assert cfg.control_cross_corner_splits == (2, 3)
        assert cfg.central_structure_splits == (4, 2)
        assert cfg.extend_splits_at_tube_boundaries is False
        assert cfg.split_tubes_in_half is True

    def test_from_yaml_true(self):
        cfg = ControlCrossSubmeshConfig.from_yaml(True)
        assert cfg is not None
        assert cfg.enabled is True
        assert cfg.control_cross_corner_splits is None
        assert cfg.extend_splits_at_tube_boundaries is True

    def test_from_yaml_false(self):
        cfg = ControlCrossSubmeshConfig.from_yaml(False)
        assert cfg is None

    def test_from_yaml_dict(self):
        raw = {
            "enabled": True,
            "control_cross_corner_splits": [2, 2],
            "central_structure_splits": [3, 2],
            "extend_splits_at_tube_boundaries": False,
            "split_tubes_in_half": True,
        }
        cfg = ControlCrossSubmeshConfig.from_yaml(raw)
        assert cfg.enabled is True
        assert cfg.control_cross_corner_splits == (2, 2)
        assert cfg.central_structure_splits == (3, 2)
        assert cfg.extend_splits_at_tube_boundaries is False
        assert cfg.split_tubes_in_half is True

    def test_from_yaml_dict_defaults(self):
        raw = {"enabled": True}
        cfg = ControlCrossSubmeshConfig.from_yaml(raw)
        assert cfg.control_cross_corner_splits is None
        assert cfg.central_structure_splits is None
        assert cfg.extend_splits_at_tube_boundaries is True
        assert cfg.split_tubes_in_half is False

    def test_from_yaml_invalid_type(self):
        with pytest.raises(TypeError, match="control_cross_submesh must be"):
            ControlCrossSubmeshConfig.from_yaml(42)

    def test_repr(self):
        cfg = ControlCrossSubmeshConfig(control_cross_corner_splits=[2, 2])
        r = repr(cfg)
        assert "ControlCrossSubmeshConfig" in r
        assert "(2, 2)" in r


# =====================================================================
#  CrossModeratorDiscretizationConfig – moderator region splits
# =====================================================================

class TestCrossModeratorDiscretizationConfig:
    def test_defaults(self):
        cdc = CrossModeratorDiscretizationConfig()
        assert cdc.narrow_gap_splits is None
        assert cdc.moderator_at_cross_corner_splits is None
        assert cdc.stub_splits is None

    def test_custom(self):
        cdc = CrossModeratorDiscretizationConfig(
            narrow_gap_splits=[4, 2],
            moderator_at_cross_corner_splits=[3, 3],
            stub_splits=[2, 1],
        )
        assert cdc.narrow_gap_splits == (4, 2)
        assert cdc.moderator_at_cross_corner_splits == (3, 3)
        assert cdc.stub_splits == (2, 1)

    def test_resolve_auto(self):
        cdc = CrossModeratorDiscretizationConfig()
        ng, mc, st = cdc.resolve(
            gap_splits=(10, 2),
            lattice_pitch=10.0,
            wide_gap_width=1.0,
            narrow_gap_width=0.5,
            cross_corner_dims=(0.5, 0.5),
            stub_dims=(5.0, 0.25),
        )
        assert isinstance(ng, tuple) and len(ng) == 2
        assert isinstance(mc, tuple) and len(mc) == 2
        assert isinstance(st, tuple) and len(st) == 2
        # All auto-computed values must be >= 1
        assert all(v >= 1 for v in ng + mc + st)

    def test_resolve_explicit(self):
        cdc = CrossModeratorDiscretizationConfig(
            narrow_gap_splits=[4, 2],
            moderator_at_cross_corner_splits=[3, 3],
            stub_splits=[2, 1],
        )
        ng, mc, st = cdc.resolve(
            gap_splits=(10, 2),
            lattice_pitch=10.0,
            wide_gap_width=1.0,
            narrow_gap_width=0.5,
            cross_corner_dims=(0.5, 0.5),
            stub_dims=(5.0, 0.25),
        )
        assert ng == (4, 2)
        assert mc == (3, 3)
        assert st == (2, 1)

    def test_repr(self):
        cdc = CrossModeratorDiscretizationConfig(
            moderator_at_cross_corner_splits=[3, 3],
        )
        r = repr(cdc)
        assert "CrossModeratorDiscretizationConfig" in r
        assert "(3, 3)" in r


# =====================================================================
#  BoxDiscretizationConfig – control_cross_submesh integration
# =====================================================================

class TestBoxDiscretizationConfigControlCrossSubmesh:
    def test_no_control_cross_submesh(self):
        bdc = BoxDiscretizationConfig()
        assert bdc.control_cross_submesh is None

    def test_control_cross_submesh_config_object(self):
        ccs = ControlCrossSubmeshConfig(control_cross_corner_splits=[2, 2])
        bdc = BoxDiscretizationConfig(control_cross_submesh=ccs)
        assert bdc.control_cross_submesh is ccs

    def test_control_cross_submesh_with_cross_moderator(self):
        cmd = CrossModeratorDiscretizationConfig(
            moderator_at_cross_corner_splits=[3, 3],
        )
        ccs = ControlCrossSubmeshConfig(
            control_cross_corner_splits=[2, 2],
        )
        bdc = BoxDiscretizationConfig(
            cross_moderator_discretization=cmd,
            control_cross_submesh=ccs,
        )
        assert bdc.cross_moderator_discretization is cmd
        assert bdc.control_cross_submesh is ccs


# =====================================================================
#  CalculationStep – construction and validation
# =====================================================================

class TestCalculationStep:

    def test_basic_ssh_step(self):
        step = CalculationStep(
            name="SSH",
            step_type="self_shielding",
            self_shielding_module="USS",
            self_shielding_method="RSE",
            spatial_method="CP",
            tracking="TISO",
        )
        assert step.name == "SSH"
        assert step.step_type == "self_shielding"
        assert step.self_shielding_module == "USS"
        assert step.self_shielding_method == "RSE"
        assert step.spatial_method == "CP"
        assert step.tracking == "TISO"
        assert step.flux_level is None
        assert step.radial_scheme == "Santamarina"
        assert step.sectorization_enabled is False
        assert step.export_macros is False

    def test_basic_flux_step(self):
        step = CalculationStep(
            name="FLUX",
            step_type="flux",
            spatial_method="IC",
            tracking="TISO",
            flux_level=1,
            export_macros=True,
        )
        assert step.step_type == "flux"
        assert step.self_shielding_module is None
        assert step.self_shielding_method is None
        assert step.flux_level == 1
        assert step.export_macros is True

    def test_invalid_step_type(self):
        with pytest.raises(ValueError, match="Invalid step_type"):
            CalculationStep(
                name="BAD", step_type="depletion",
                self_shielding_module="USS", self_shielding_method="RSE",
                spatial_method="CP",
            )

    def test_invalid_self_shielding_module(self):
        with pytest.raises(ValueError, match="Invalid self_shielding_module"):
            CalculationStep(
                name="BAD", step_type="self_shielding",
                self_shielding_module="MAGIC", self_shielding_method="RSE",
                spatial_method="CP",
            )

    def test_invalid_self_shielding_method(self):
        with pytest.raises(ValueError, match="Invalid self_shielding_method"):
            CalculationStep(
                name="BAD", step_type="self_shielding",
                self_shielding_module="USS", self_shielding_method="MAGIC",
                spatial_method="CP",
            )

    def test_self_shielding_module_required_for_ssh_step(self):
        with pytest.raises(ValueError, match="self_shielding_module is required"):
            CalculationStep(
                name="BAD", step_type="self_shielding",
                self_shielding_method="RSE",
                spatial_method="CP",
            )

    def test_self_shielding_method_required_for_ssh_step(self):
        with pytest.raises(ValueError, match="self_shielding_method is required"):
            CalculationStep(
                name="BAD", step_type="self_shielding",
                self_shielding_module="USS",
                spatial_method="CP",
            )

    def test_self_shielding_module_not_allowed_for_flux_step(self):
        with pytest.raises(ValueError, match="self_shielding_module should not be set for flux"):
            CalculationStep(
                name="BAD", step_type="flux",
                self_shielding_module="USS",
                spatial_method="CP",
            )

    def test_self_shielding_method_not_allowed_for_flux_step(self):
        with pytest.raises(ValueError, match="self_shielding_method should not be set for flux"):
            CalculationStep(
                name="BAD", step_type="flux",
                self_shielding_method="RSE",
                spatial_method="CP",
            )

    def test_invalid_spatial_method(self):
        with pytest.raises(ValueError, match="Invalid spatial_method"):
            CalculationStep(
                name="BAD", step_type="flux",
                spatial_method="FEM",
            )

    def test_moc_only_for_flux(self):
        with pytest.raises(ValueError, match="MOC spatial method is only available for flux"):
            CalculationStep(
                name="BAD", step_type="self_shielding",
                self_shielding_module="USS", self_shielding_method="RSE",
                spatial_method="MOC",
            )

    def test_ic_requires_tiso(self):
        with pytest.raises(ValueError, match="only supports 'TISO'"):
            CalculationStep(
                name="BAD", step_type="flux",
                spatial_method="IC",
                tracking="TSPC",
            )

    def test_invalid_tracking(self):
        with pytest.raises(ValueError, match="Invalid tracking"):
            CalculationStep(
                name="BAD", step_type="flux",
                spatial_method="CP",
                tracking="ANISO",
            )

    def test_invalid_radial_scheme(self):
        with pytest.raises(ValueError, match="Invalid radial_scheme"):
            CalculationStep(
                name="BAD", step_type="flux",
                spatial_method="CP",
                radial_scheme="random",
            )

    def test_moc_with_tspc_valid(self):
        step = CalculationStep(
            name="FLUX_MOC", step_type="flux",
            spatial_method="MOC",
            tracking="TSPC",
            polar_angles_quadrature="GAUS",
            number_of_polar_angles=4,
        )
        assert step.spatial_method == "MOC"
        assert step.tracking == "TSPC"
        assert step.polar_angles_quadrature == "GAUS"
        assert step.number_of_polar_angles == 4

    def test_sectorization_query_disabled(self):
        step = CalculationStep(
            name="S", step_type="flux",
            spatial_method="CP",
            sectorization_enabled=False,
            fuel_sectors=SectorConfig(sectors=[8]),
        )
        assert step.get_sectorization_for_pin("ROD1") is None
        assert step.get_water_rod_sectorization() is None

    def test_sectorization_query_enabled(self):
        fuel_sc = SectorConfig(sectors=[1, 1, 8])
        gd_sc = SectorConfig(sectors=[1, 1, 1, 8])
        wr_sc = SectorConfig(sectors=[1, 8])
        step = CalculationStep(
            name="S", step_type="flux",
            spatial_method="IC",
            sectorization_enabled=True,
            fuel_sectors=fuel_sc,
            gd_sectors=gd_sc,
            water_rod_sectors=wr_sc,
        )
        assert step.get_sectorization_for_pin("ROD1") is fuel_sc
        assert step.get_sectorization_for_pin("ROD5G", isGd=True) is gd_sc
        assert step.get_water_rod_sectorization() is wr_sc

    def test_repr(self):
        step = CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
        )
        r = repr(step)
        assert "SSH" in r
        assert "self_shielding" in r

    def test_radial_overrides(self):
        step = CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
            radial_scheme="Santamarina",
            radial_overrides={
                "Gd": {"scheme": "automatic", "params": {"num_radial_zones": 6}},
                "ROD1": {"scheme": "user_defined", "params": {"user_defined_radii": [0.1, 0.3, 0.438]}},
            },
        )
        # Create a mock pin to test resolution
        class MockPin:
            def __init__(self, rod_ID, isGd):
                self.rod_ID = rod_ID
                self.isGd = isGd

        # Specific rod override
        scheme, params = step._resolve_radial_config(MockPin("ROD1", False))
        assert scheme == "user_defined"
        assert params["user_defined_radii"] == [0.1, 0.3, 0.438]

        # Gd category override
        scheme, params = step._resolve_radial_config(MockPin("ROD5G", True))
        assert scheme == "automatic"
        assert params["num_radial_zones"] == 6

        # Fallback to step default
        scheme, params = step._resolve_radial_config(MockPin("ROD2", False))
        assert scheme == "Santamarina"


# =====================================================================
#  DragonCalculationScheme – programmatic construction
# =====================================================================

class TestDragonCalculationScheme:

    def test_empty_scheme(self):
        scheme = DragonCalculationScheme(name="empty")
        assert scheme.name == "empty"
        assert len(scheme.steps) == 0
        assert scheme.get_self_shielding_steps() == []
        assert scheme.get_flux_steps() == []

    def test_add_step(self):
        scheme = DragonCalculationScheme()
        step = CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
        )
        scheme.add_step(step)
        assert len(scheme.steps) == 1
        assert scheme.steps[0] is step

    def test_add_step_type_check(self):
        scheme = DragonCalculationScheme()
        with pytest.raises(TypeError, match="Expected CalculationStep"):
            scheme.add_step("not_a_step")

    def test_get_step_by_name(self):
        scheme = DragonCalculationScheme()
        step = CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
        )
        scheme.add_step(step)
        assert scheme.get_step("SSH") is step

    def test_get_step_not_found(self):
        scheme = DragonCalculationScheme()
        with pytest.raises(KeyError, match="No step named"):
            scheme.get_step("NONEXISTENT")

    def test_get_steps_by_type(self):
        scheme = DragonCalculationScheme()
        ssh = CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
        )
        flux1 = CalculationStep(
            name="FLUX_L1", step_type="flux",
            spatial_method="IC",
            flux_level=1,
        )
        flux2 = CalculationStep(
            name="FLUX_L2", step_type="flux",
            spatial_method="MOC",
            tracking="TSPC", flux_level=2,
            polar_angles_quadrature="GAUS",
            number_of_polar_angles=4,
        )
        scheme.add_step(ssh)
        scheme.add_step(flux2)  # Add L2 first
        scheme.add_step(flux1)  # Add L1 second

        ssh_steps = scheme.get_self_shielding_steps()
        assert len(ssh_steps) == 1
        assert ssh_steps[0].name == "SSH"

        flux_steps = scheme.get_flux_steps()
        assert len(flux_steps) == 2
        # Should be sorted by flux_level
        assert flux_steps[0].flux_level == 1
        assert flux_steps[1].flux_level == 2

    def test_repr(self):
        scheme = DragonCalculationScheme(name="test_scheme")
        scheme.add_step(CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
        ))
        r = repr(scheme)
        assert "test_scheme" in r
        assert "SSH" in r

    def test_summary(self):
        scheme = DragonCalculationScheme.preset("BWR_fine_1L")
        summary = scheme.summary()
        assert "SSH" in summary
        assert "FLUX" in summary
        assert "self_shielding" in summary
        assert "Santamarina" in summary


# =====================================================================
#  DragonCalculationScheme – presets
# =====================================================================

class TestPresets:

    def test_BWR_fine_1L(self):
        scheme = DragonCalculationScheme.preset("BWR_fine_1L")
        assert scheme.name == "BWR_fine_1L"
        assert len(scheme.steps) == 2

        ssh = scheme.get_step("SSH")
        assert ssh.step_type == "self_shielding"
        assert ssh.self_shielding_module == "USS"
        assert ssh.self_shielding_method == "RSE"
        assert ssh.spatial_method == "IC"
        assert ssh.tracking == "TISO"
        assert ssh.radial_scheme == "Santamarina"
        assert ssh.sectorization_enabled is False
        assert ssh.export_macros is True

        flux = scheme.get_step("FLUX")
        assert flux.step_type == "flux"
        assert flux.self_shielding_module is None
        assert flux.self_shielding_method is None
        assert flux.spatial_method == "MOC"
        assert flux.tracking == "TSPC"
        assert flux.sectorization_enabled is True
        assert flux.export_macros is False
        assert flux.fuel_sectors is not None
        assert flux.fuel_sectors.windmill is True
        assert flux.gd_sectors is not None
        assert flux.water_rod_sectors is not None
        assert flux.box_discretization is not None
        assert flux.box_discretization.enabled is True

    def test_BWR_2L(self):
        scheme = DragonCalculationScheme.preset("BWR_2L")
        assert scheme.name == "BWR_2L"
        assert len(scheme.steps) == 3

        ssh = scheme.get_step("SSH")
        assert ssh.step_type == "self_shielding"
        assert ssh.self_shielding_module == "USS"
        assert ssh.self_shielding_method == "RSE"

        flux_steps = scheme.get_flux_steps()
        assert len(flux_steps) == 2
        assert flux_steps[0].flux_level == 1
        assert flux_steps[0].spatial_method == "IC"
        assert flux_steps[0].self_shielding_module is None
        assert flux_steps[1].flux_level == 2
        assert flux_steps[1].spatial_method == "MOC"
        assert flux_steps[1].tracking == "TSPC"
        assert flux_steps[1].box_discretization is not None
        assert flux_steps[1].box_discretization.enabled is True

    def test_BWR_CP(self):
        scheme = DragonCalculationScheme.preset("BWR_CP")
        assert scheme.name == "BWR_CP"
        assert len(scheme.steps) == 2

        ssh = scheme.get_step("SSH")
        assert ssh.step_type == "self_shielding"
        assert ssh.self_shielding_module == "USS"
        assert ssh.self_shielding_method == "PT"
        assert ssh.spatial_method == "CP"
        assert ssh.sectorization_enabled is False

        flux = scheme.get_step("FLUX")
        assert flux.step_type == "flux"
        assert flux.self_shielding_module is None
        assert flux.spatial_method == "CP"
        assert flux.tracking == "TSPC"
        assert flux.sectorization_enabled is False

    def test_invalid_preset(self):
        with pytest.raises(ValueError, match="Unknown preset"):
            DragonCalculationScheme.preset("nonexistent")


# =====================================================================
#  DragonCalculationScheme – YAML construction
# =====================================================================

class TestFromYAML:

    def test_from_yaml_GE14(self):
        scheme = DragonCalculationScheme.from_yaml(_CALC_SCHEME_YAML)
        assert scheme.name == "GE14_standard"
        assert len(scheme.steps) == 2

        ssh = scheme.get_step("SSH")
        assert ssh.step_type == "self_shielding"
        assert ssh.self_shielding_module == "USS"
        assert ssh.self_shielding_method == "RSE"
        assert ssh.spatial_method == "CP"
        assert ssh.sectorization_enabled is False

        flux = scheme.get_step("FLUX")
        assert flux.step_type == "flux"
        assert flux.self_shielding_module is None
        assert flux.self_shielding_method is None
        assert flux.spatial_method == "IC"
        assert flux.tracking == "TISO"
        assert flux.export_macros is True
        assert flux.sectorization_enabled is True
        assert flux.fuel_sectors.windmill is True
        assert flux.fuel_sectors.sectors == [1, 1, 1, 1, 1, 1, 8]
        assert flux.gd_sectors.sectors == [1, 1, 1, 1, 1, 1, 1, 1, 8]
        assert flux.water_rod_sectors.sectors == [1, 1, 8]

    def test_from_yaml_roundtrip(self):
        """Write a scheme to YAML, read it back, and verify consistency."""
        scheme_data = {
            "DRAGON_CALCULATION_SCHEME": {
                "name": "roundtrip_test",
                "steps": [
                    {
                        "name": "SSH",
                        "step_type": "self_shielding",
                        "self_shielding_module": "USS",
                        "self_shielding_method": "RSE",
                        "spatial_method": "CP",
                        "tracking": "TISO",
                        "radial_scheme": "Santamarina",
                        "sectorization": {"enabled": False},
                    },
                    {
                        "name": "FLUX_MOC",
                        "step_type": "flux",
                        "spatial_method": "MOC",
                        "tracking": "TSPC",
                        "flux_level": 2,
                        "radial_scheme": "automatic",
                        "radial_params": {"num_radial_zones": 3},
                        "export_macros": True,
                        "num_angles_2d": 24,
                        "line_density": 40.0,
                        "anisotropy_level": 4,
                        "polar_angles_quadrature": "GAUS",
                        "number_of_polar_angles": 4,
                        "sectorization": {
                            "enabled": True,
                            "windmill": True,
                            "fuel_pins": {
                                "sectors": [1, 1, 1, 12],
                                "angles": [0, 0, 0, 15.0],
                            },
                        },
                        "box_discretization": {
                            "enabled": True,
                            "corner_splits": [3, 3],
                            "gap_splits": [10, 1],
                        },
                    },
                ],
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(scheme_data, f)
            tmp_path = f.name

        try:
            scheme = DragonCalculationScheme.from_yaml(tmp_path)
            assert scheme.name == "roundtrip_test"
            assert len(scheme.steps) == 2

            moc_step = scheme.get_step("FLUX_MOC")
            assert moc_step.spatial_method == "MOC"
            assert moc_step.tracking == "TSPC"
            assert moc_step.flux_level == 2
            assert moc_step.radial_scheme == "automatic"
            assert moc_step.radial_params == {"num_radial_zones": 3}
            assert moc_step.sectorization_enabled is True
            assert moc_step.fuel_sectors.windmill is True
            assert moc_step.fuel_sectors.sectors == [1, 1, 1, 12]
            assert moc_step.export_macros is True
            assert moc_step.box_discretization is not None
            assert moc_step.box_discretization.enabled is True
            assert moc_step.box_discretization.corner_splits == (3, 3)
            assert moc_step.box_discretization.gap_splits == (10, 1)
        finally:
            os.unlink(tmp_path)

    def test_from_yaml_water_rod_splits(self):
        """Verify water_rods 'splits' is parsed into SectorConfig.splits."""
        scheme_data = {
            "DRAGON_CALCULATION_SCHEME": {
                "name": "wr_splits_test",
                "steps": [
                    {
                        "name": "FLUX",
                        "step_type": "flux",
                        "spatial_method": "IC",
                        "tracking": "TISO",
                        "radial_scheme": "Santamarina",
                        "export_macros": True,
                        "sectorization": {
                            "enabled": True,
                            "windmill": False,
                            "fuel_pins": {
                                "sectors": [1, 1, 1, 1, 1, 1, 8],
                                "angles": [0, 0, 0, 0, 0, 0, 22.5],
                            },
                            "water_rods": {
                                "splits": [3, 3],
                            },
                        },
                    },
                ],
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(scheme_data, f)
            tmp_path = f.name

        try:
            scheme = DragonCalculationScheme.from_yaml(tmp_path)
            flux = scheme.get_step("FLUX")
            assert flux.water_rod_sectors is not None
            assert flux.water_rod_sectors.splits == (3, 3)
            # sectors/angles should be empty since only splits was specified
            assert flux.water_rod_sectors.sectors == []
            assert flux.water_rod_sectors.angles == []
        finally:
            os.unlink(tmp_path)

    def test_from_yaml_water_rod_splits_and_sectors_warns(self):
        """Warn when water_rods block has both splits and sectors/angles."""
        import warnings
        scheme_data = {
            "DRAGON_CALCULATION_SCHEME": {
                "name": "wr_mixed_test",
                "steps": [
                    {
                        "name": "FLUX",
                        "step_type": "flux",
                        "spatial_method": "IC",
                        "tracking": "TISO",
                        "radial_scheme": "Santamarina",
                        "sectorization": {
                            "enabled": True,
                            "water_rods": {
                                "splits": [3, 3],
                                "sectors": [1, 1, 8],
                                "angles": [0, 0, 22.5],
                            },
                        },
                    },
                ],
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(scheme_data, f)
            tmp_path = f.name

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                scheme = DragonCalculationScheme.from_yaml(tmp_path)
                wr_warnings = [x for x in w if "splits" in str(x.message)
                               and "sectors" in str(x.message)]
                assert len(wr_warnings) >= 1, \
                    "Expected a warning about both splits and sectors being set"
            # Both should still be stored
            flux = scheme.get_step("FLUX")
            assert flux.water_rod_sectors.splits == (3, 3)
            assert flux.water_rod_sectors.sectors == [1, 1, 8]
        finally:
            os.unlink(tmp_path)

    def test_from_yaml_additional_radial_splits_int(self):
        """Verify additional_radial_splits_in_moderator (int) is parsed
        into SectorConfig."""
        scheme_data = {
            "DRAGON_CALCULATION_SCHEME": {
                "name": "wr_radial_int_test",
                "steps": [
                    {
                        "name": "FLUX",
                        "step_type": "flux",
                        "spatial_method": "MOC",
                        "tracking": "TSPC",
                        "radial_scheme": "Santamarina",
                        "polar_angles_quadrature": "GAUS",
                        "number_of_polar_angles": 4,
                        "sectorization": {
                            "enabled": True,
                            "windmill": True,
                            "fuel_pins": {
                                "sectors": [8, 8, 8, 8, 8, 8, 8],
                                "angles": [22.5] * 7,
                            },
                            "water_rods": {
                                "additional_radial_splits_in_moderator": 3,
                                "sectors": [16, 16, 16],
                                "angles": [0, 0, 0],
                            },
                        },
                    },
                ],
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(scheme_data, f)
            tmp_path = f.name

        try:
            scheme = DragonCalculationScheme.from_yaml(tmp_path)
            flux = scheme.get_step("FLUX")
            wr = flux.water_rod_sectors
            assert wr is not None
            assert wr.additional_radial_splits_in_moderator == 3
            assert wr.sectors == [16, 16, 16]
            # Verify expanded sectors: 2 extra + 3 base = 5
            expanded_s, expanded_a = wr.expanded_sectors_and_angles(1.0)
            assert len(expanded_s) == 5
            assert expanded_s == [16, 16, 16, 16, 16]
        finally:
            os.unlink(tmp_path)

    def test_from_yaml_additional_radial_splits_list(self):
        """Verify additional_radial_splits_in_moderator (list) is parsed
        into SectorConfig."""
        scheme_data = {
            "DRAGON_CALCULATION_SCHEME": {
                "name": "wr_radial_list_test",
                "steps": [
                    {
                        "name": "FLUX",
                        "step_type": "flux",
                        "spatial_method": "MOC",
                        "tracking": "TSPC",
                        "radial_scheme": "Santamarina",
                        "polar_angles_quadrature": "GAUS",
                        "number_of_polar_angles": 4,
                        "sectorization": {
                            "enabled": True,
                            "water_rods": {
                                "additional_radial_splits_in_moderator": [0.3, 0.6],
                                "sectors": [8, 8, 8],
                                "angles": [0, 0, 0],
                            },
                        },
                    },
                ],
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(scheme_data, f)
            tmp_path = f.name

        try:
            scheme = DragonCalculationScheme.from_yaml(tmp_path)
            flux = scheme.get_step("FLUX")
            wr = flux.water_rod_sectors
            assert wr.additional_radial_splits_in_moderator == [0.3, 0.6]
            radii = wr.resolve_water_rod_radii(1.0)
            assert radii == [0.3, 0.6]
            # Expanded: 2 extra + 3 base = 5
            expanded_s, _ = wr.expanded_sectors_and_angles(1.0)
            assert len(expanded_s) == 5
        finally:
            os.unlink(tmp_path)

    def test_from_yaml_additional_radial_splits_default(self):
        """When additional_radial_splits_in_moderator is absent, default is 1."""
        scheme_data = {
            "DRAGON_CALCULATION_SCHEME": {
                "name": "wr_radial_default_test",
                "steps": [
                    {
                        "name": "FLUX",
                        "step_type": "flux",
                        "spatial_method": "IC",
                        "tracking": "TISO",
                        "radial_scheme": "Santamarina",
                        "sectorization": {
                            "enabled": True,
                            "water_rods": {
                                "sectors": [1, 1, 8],
                                "angles": [0, 0, 22.5],
                            },
                        },
                    },
                ],
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(scheme_data, f)
            tmp_path = f.name

        try:
            scheme = DragonCalculationScheme.from_yaml(tmp_path)
            flux = scheme.get_step("FLUX")
            wr = flux.water_rod_sectors
            assert wr.additional_radial_splits_in_moderator == None
            assert wr.resolve_water_rod_radii(1.0) == []
            s, a = wr.expanded_sectors_and_angles(1.0)
            assert s == [1, 1, 8]
        finally:
            os.unlink(tmp_path)

    def test_from_yaml_radial_overrides(self):
        """Check that per-rod radial overrides are parsed from YAML."""
        scheme_data = {
            "DRAGON_CALCULATION_SCHEME": {
                "name": "override_test",
                "steps": [
                    {
                        "name": "SSH",
                        "step_type": "self_shielding",
                        "self_shielding_module": "USS",
                        "self_shielding_method": "RSE",
                        "spatial_method": "CP",
                        "radial_scheme": "Santamarina",
                        "radial_overrides": {
                            "Gd": {
                                "scheme": "automatic",
                                "params": {"num_radial_zones": 6},
                            },
                        },
                        "sectorization": {"enabled": False},
                    },
                ],
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(scheme_data, f)
            tmp_path = f.name

        try:
            scheme = DragonCalculationScheme.from_yaml(tmp_path)
            ssh = scheme.get_step("SSH")
            assert "Gd" in ssh.radial_overrides
            assert ssh.radial_overrides["Gd"]["scheme"] == "automatic"
            assert ssh.radial_overrides["Gd"]["params"]["num_radial_zones"] == 6
        finally:
            os.unlink(tmp_path)


# =====================================================================
#  CalculationStep.apply_radii – integration with assembly model
# =====================================================================

class TestApplyRadii:

    def test_apply_santamarina_radii(self):
        """apply_radii with Santamarina scheme should produce 4 fuel zones
        for UOX and 6 for Gd pins on the simplified assembly."""
        assembly = _build_simple_assembly()

        step = CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
            radial_scheme="Santamarina",
        )
        step.apply_radii(assembly)

        for row in assembly.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    pin.count_number_of_fuel_radial_self_shielding_zones()
                    if pin.isGd:
                        assert pin.number_of_self_shielding_fuel_zones == 6, \
                            f"Gd pin should have 6 Santamarina zones, got {pin.number_of_self_shielding_fuel_zones}"
                    else:
                        assert pin.number_of_self_shielding_fuel_zones == 4, \
                            f"UOX pin should have 4 Santamarina zones, got {pin.number_of_self_shielding_fuel_zones}"

    def test_apply_automatic_radii(self):
        """apply_radii with automatic scheme and 3 zones."""
        assembly = _build_simple_assembly()

        step = CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
            radial_scheme="automatic",
            radial_params={"num_radial_zones": 3},
        )
        step.apply_radii(assembly)

        for row in assembly.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    pin.count_number_of_fuel_radial_self_shielding_zones()
                    assert pin.number_of_self_shielding_fuel_zones == 3, \
                        f"Expected 3 automatic zones, got {pin.number_of_self_shielding_fuel_zones}"

    def test_apply_radii_with_gd_override(self):
        """Gd override: Santamarina for fuel, automatic(6) for Gd."""
        assembly = _build_simple_assembly()

        step = CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
            radial_scheme="Santamarina",
            radial_overrides={
                "Gd": {"scheme": "automatic", "params": {"num_radial_zones": 6}},
            },
        )
        step.apply_radii(assembly)

        for row in assembly.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    pin.count_number_of_fuel_radial_self_shielding_zones()
                    if pin.isGd:
                        assert pin.number_of_self_shielding_fuel_zones == 6
                    else:
                        assert pin.number_of_self_shielding_fuel_zones == 4

    def test_apply_radii_requires_lattice(self):
        """apply_radii should fail if lattice hasn't been built."""
        assembly = CartesianAssemblyModel(
            name="no_lattice", tdt_file="dummy.tdt",
            geometry_description_yaml=_GEOMETRY_SIMPLE_YAML,
        )
        step = CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
        )
        with pytest.raises(RuntimeError, match="lattice has not been built"):
            step.apply_radii(assembly)


# =====================================================================
#  YAML key validation
# =====================================================================

class TestYAMLValidation:

    def test_wrong_key_detected(self):
        """The YAML validator should reject 'water_rod_ids' (wrong key)."""
        geometry_data = {
            "PIN_GEOMETRY": {
                "fuel_radius": 0.438,
                "pin_pitch": 1.3,
            },
            "ASSEMBLY_GEOMETRY": {
                "lattice_type": "Cartesian",
                "lattice_description": [["ROD1"]],
                "assembly_pitch": 1.3,
                "water_rod_ids": ["WROD"],  # Wrong key!
            },
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(geometry_data, f)
            tmp_path = f.name

        try:
            with pytest.raises(ValueError, match="Unrecognised key.*water_rod_ids"):
                CartesianAssemblyModel(
                    name="bad_yaml", tdt_file="dummy.tdt",
                    geometry_description_yaml=tmp_path,
                )
        finally:
            os.unlink(tmp_path)

    def test_correct_key_accepted(self):
        """The YAML validator should accept 'non_fuel_rod_ids'."""
        geometry_data = {
            "PIN_GEOMETRY": {
                "fuel_radius": 0.438,
                "pin_pitch": 1.3,
            },
            "ASSEMBLY_GEOMETRY": {
                "lattice_type": "Cartesian",
                "lattice_description": [["ROD1"]],
                "assembly_pitch": 1.3,
                "non_fuel_rod_ids": ["WROD"],
            },
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(geometry_data, f)
            tmp_path = f.name

        try:
            # Should not raise
            assembly = CartesianAssemblyModel(
                name="good_yaml", tdt_file="dummy.tdt",
                geometry_description_yaml=tmp_path,
            )
            assert assembly.non_fuel_rod_ids == ["WROD"]
        finally:
            os.unlink(tmp_path)

    def test_missing_assembly_geometry_section(self):
        """Should raise if ASSEMBLY_GEOMETRY section is missing entirely."""
        geometry_data = {
            "PIN_GEOMETRY": {
                "fuel_radius": 0.438,
                "pin_pitch": 1.3,
            },
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(geometry_data, f)
            tmp_path = f.name

        try:
            with pytest.raises(ValueError, match="missing.*ASSEMBLY_GEOMETRY"):
                CartesianAssemblyModel(
                    name="no_asm", tdt_file="dummy.tdt",
                    geometry_description_yaml=tmp_path,
                )
        finally:
            os.unlink(tmp_path)

    def test_missing_required_key(self):
        """Should raise if a required key (lattice_description) is missing."""
        geometry_data = {
            "ASSEMBLY_GEOMETRY": {
                "assembly_pitch": 1.3,
                # lattice_description is missing
            },
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(geometry_data, f)
            tmp_path = f.name

        try:
            with pytest.raises(ValueError, match="Missing required key"):
                CartesianAssemblyModel(
                    name="no_lattice_desc", tdt_file="dummy.tdt",
                    geometry_description_yaml=tmp_path,
                )
        finally:
            os.unlink(tmp_path)
