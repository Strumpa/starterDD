# Tests for thermal scattering integration
# Date : 26/02/2026
# R.Guasch
#
# Covers:
#   - Composition.setTherm (auto-detect + manual dict)
#   - DEFAULT_THERMAL_SCATTERING registry
#   - get_xs_suffix (closest-inferior logic)
#   - get_therm_suffix / get_therm_interpolation_suffixes
#   - S2_ThermalScattering (single library + interpolation)
#   - S2_Settings.set_nuclear_data_evaluation
#   - Evaluation-dependent suffix maps

import pytest

from starterDD.MaterialProperties.material_mixture import (
    Composition,
    DEFAULT_THERMAL_SCATTERING,
)

from starterDD.InterfaceToDD.serpent2_cards import (
    get_xs_suffix,
    get_therm_suffix,
    get_therm_interpolation_suffixes,
    S2_ThermalScattering,
    S2_Settings,
    EVALUATION_XS_SUFFIX_MAPS,
    EVALUATION_THERM_SUFFIX_MAPS,
    TEMPERATURE_TO_XS_SUFFIX,
)


# ---------------------------------------------------------------------------
#  DEFAULT_THERMAL_SCATTERING registry
# ---------------------------------------------------------------------------
class TestDefaultThermalScatteringRegistry:
    def test_H1_entry_exists(self):
        assert "H1" in DEFAULT_THERMAL_SCATTERING
        entry = DEFAULT_THERMAL_SCATTERING["H1"]
        assert entry["dragon_alias"] == "H1_H2O"
        assert entry["serpent2_therm_name"] == "lwtr"
        assert entry["serpent2_zaid"] == "1001"

    def test_H2_entry_exists(self):
        assert "H2" in DEFAULT_THERMAL_SCATTERING
        entry = DEFAULT_THERMAL_SCATTERING["H2"]
        assert entry["dragon_alias"] == "D2_D2O"
        assert entry["serpent2_therm_name"] == "hwtr"
        assert entry["serpent2_zaid"] == "1002"


# ---------------------------------------------------------------------------
#  Composition.setTherm
# ---------------------------------------------------------------------------
class TestCompositionSetTherm:
    def test_auto_detect_H1_in_water(self):
        """therm=True should auto-detect H1 from DEFAULT_THERMAL_SCATTERING."""
        comp = Composition("water", {"1001": 4.9e-2, "8016": 2.5e-2})
        comp.setTherm(True)
        assert comp.therm is True
        assert len(comp.therm_data) == 1
        assert comp.therm_data[0]["isotope"] == "H1"
        assert comp.therm_data[0]["dragon_alias"] == "H1_H2O"

    def test_auto_detect_no_match(self):
        """Composition without matching isotopes → therm=True but empty data."""
        comp = Composition("zircaloy", {"40090": 2.2e-2, "40091": 4.8e-3})
        comp.setTherm(True)
        assert comp.therm is True
        assert comp.therm_data == []

    def test_manual_dict(self):
        """Explicit dict overrides the default registry."""
        comp = Composition("graphite", {"6012": 8.0e-2})
        comp.setTherm({
            "C12": {
                "dragon_alias": "C12_GRAPH",
                "serpent2_therm_name": "grph",
                "serpent2_zaid": "6012",
            }
        })
        assert comp.therm is True
        assert len(comp.therm_data) == 1
        assert comp.therm_data[0]["isotope"] == "C12"
        assert comp.therm_data[0]["serpent2_therm_name"] == "grph"

    def test_setTherm_false_clears(self):
        comp = Composition("water", {"1001": 4.9e-2, "8016": 2.5e-2})
        comp.setTherm(True)
        assert comp.therm is True
        comp.setTherm(False)
        assert comp.therm is False
        assert comp.therm_data == []


# ---------------------------------------------------------------------------
#  get_xs_suffix — closest-inferior logic
# ---------------------------------------------------------------------------
class TestGetXsSuffix:
    def test_exact_match(self):
        assert get_xs_suffix(900.0) == '.09c'

    def test_inferior_temperature(self):
        """559 K is between 550 and 900 → should pick 550 (inferior)."""
        assert get_xs_suffix(559.0) == '.05c'

    def test_above_highest(self):
        """Temperature above all tabulated → should pick highest inferior."""
        assert get_xs_suffix(2000.0) == '.20c'

    def test_below_lowest_warns(self):
        """Temperature below all tabulated → warning + lowest available."""
        with pytest.warns(match="below all tabulated"):
            suffix = get_xs_suffix(100.0)
        assert suffix == '.02c'

    def test_custom_suffix_map(self):
        custom = {300.0: '.30c', 600.0: '.60c', 900.0: '.90c'}
        assert get_xs_suffix(750.0, suffix_map=custom) == '.60c'


# ---------------------------------------------------------------------------
#  get_therm_suffix / get_therm_interpolation_suffixes
# ---------------------------------------------------------------------------
class TestThermSuffix:
    def test_exact_temperature(self):
        lo, hi = get_therm_interpolation_suffixes(294.0, evaluation="endfb8r1")
        assert hi is None  # exact match
        assert lo == '.00t'

    def test_interpolation_between_points(self):
        """559 K is between 500 (.04t) and 600 (.05t)."""
        lo, hi = get_therm_interpolation_suffixes(559.0, evaluation="endfb8r1")
        assert lo == '.05t'
        assert hi == '.06t'

    def test_above_all_points(self):
        """Above all tabulated → lo is the highest, hi is None."""
        lo, hi = get_therm_interpolation_suffixes(1200.0, evaluation="endfb8r1")
        assert hi is None

    def test_get_therm_suffix_returns_inferior(self):
        suffix = get_therm_suffix(559.0, evaluation="endfb8r1")
        assert suffix == '.05t'


# ---------------------------------------------------------------------------
#  S2_ThermalScattering
# ---------------------------------------------------------------------------
class TestS2ThermalScattering:
    def test_single_library_card(self):
        ts = S2_ThermalScattering(name="lwtr", library_name="lwtr.06t")
        assert ts.format_card() == "therm lwtr  lwtr.06t"

    def test_interpolation_card(self):
        ts = S2_ThermalScattering(
            name="lwtr",
            library_name=None,
            interpolation_temperature=559.0,
            lo_library="lwtr.06t",
            hi_library="lwtr.07t",
        )
        expected = "therm lwtr 559.0 lwtr.06t lwtr.07t"
        assert ts.format_card() == expected

    def test_from_temperature_exact(self):
        ts = S2_ThermalScattering.from_temperature(
            "lwtr", 294.0, evaluation="endfb8r1"
        )
        assert ts.library_name == "lwtr.00t"
        assert ts.interpolation_temperature is None

    def test_from_temperature_interpolation(self):
        ts = S2_ThermalScattering.from_temperature(
            "lwtr", 559.0, evaluation="endfb8r1"
        )
        assert ts.interpolation_temperature == 559.0
        assert ts.lo_library == "lwtr.05t"
        assert ts.hi_library == "lwtr.06t"


# ---------------------------------------------------------------------------
#  S2_Settings — nuclear data evaluation
# ---------------------------------------------------------------------------
class TestS2SettingsEvaluation:
    def test_set_valid_evaluation(self):
        s = S2_Settings()
        s.set_nuclear_data_evaluation("endfb8r1")
        assert s.nuclear_data_evaluation == "endfb8r1"

    def test_get_xs_suffix_map_with_evaluation(self):
        s = S2_Settings()
        s.set_nuclear_data_evaluation("jeff311")
        suffix_map = s.get_xs_suffix_map()
        assert suffix_map is EVALUATION_XS_SUFFIX_MAPS["jeff311"]

    def test_get_xs_suffix_map_without_evaluation(self):
        s = S2_Settings()
        suffix_map = s.get_xs_suffix_map()
        assert suffix_map is TEMPERATURE_TO_XS_SUFFIX

    def test_unknown_evaluation_warns(self):
        s = S2_Settings()
        with pytest.warns(match="not found"):
            s.set_nuclear_data_evaluation("unknown_lib")
        assert s.nuclear_data_evaluation == "unknown_lib"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
