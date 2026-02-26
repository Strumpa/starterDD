# Tests for isotopic density conversion utilities
# Date : 26/02/2026
# R.Guasch
#
# Covers:
#   - fractions_to_iso_densities (atomic_density & mass_density modes)
#   - get_isotope_atomic_mass
#   - parse_all_compositions_from_yaml (round-trip with the GE-14 YAML)

import math
import pytest

from conftest import GE14_COMPOSITIONS_YAML

from starterDD.MaterialProperties.material_mixture import (
    fractions_to_iso_densities,
    get_isotope_atomic_mass,
    parse_all_compositions_from_yaml,
    AVOGADRO,
    CM2_TO_BARN,
    DEFAULT_THERMAL_SCATTERING,
)


# ---------------------------------------------------------------------------
#  get_isotope_atomic_mass
# ---------------------------------------------------------------------------
class TestGetIsotopeAtomicMass:
    def test_hydrogen(self):
        assert get_isotope_atomic_mass("1001") == 1.0

    def test_uranium235(self):
        assert get_isotope_atomic_mass("92235") == 235.0

    def test_oxygen16(self):
        assert get_isotope_atomic_mass("8016") == 16.0

    def test_natural_element_raises(self):
        """A=0 (natural element) should raise since individual isotopes are required."""
        with pytest.raises(ValueError, match="Natural-element"):
            get_isotope_atomic_mass("92000")


# ---------------------------------------------------------------------------
#  fractions_to_iso_densities — atomic_density mode
# ---------------------------------------------------------------------------
class TestAtomicDensityMode:
    """N_i = f_i * N_tot"""

    def test_water_atomic_fractions(self):
        """Reproduce the GE-14 moderator entry:
        atomic_density = 7.38974602E-2 atoms/barn·cm
        composition: H-1 = 2/3, O-16 = 1/3
        """
        N_tot = 7.38974602e-2
        comp = {"1001": 2.0 / 3.0, "8016": 1.0 / 3.0}
        result = fractions_to_iso_densities(comp, "atomic_density", N_tot)

        assert result["1001"] == pytest.approx(N_tot * 2.0 / 3.0, rel=1e-10)
        assert result["8016"] == pytest.approx(N_tot * 1.0 / 3.0, rel=1e-10)

    def test_fractions_must_sum_to_one(self):
        with pytest.raises(ValueError, match="fractions sum"):
            fractions_to_iso_densities({"1001": 0.5}, "atomic_density", 1.0)


# ---------------------------------------------------------------------------
#  fractions_to_iso_densities — mass_density mode
# ---------------------------------------------------------------------------
class TestMassDensityMode:
    """N_i = rho * N_A * w_i / (A_i * 1e24)"""

    def test_water_mass_fractions(self):
        """Water at ~1 g/cm³ with mass fractions of H and O."""
        rho = 1.0  # g/cm³
        M_H = 1.00794
        M_O = 15.9994
        M_H2O = 2.0 * M_H + M_O  # ~18.01528
        w_H = (2.0 * M_H) / M_H2O
        w_O = M_O / M_H2O

        # Using A (mass number) as atomic mass (our current approximation)
        A_H = 1.0
        A_O = 16.0

        comp = {"1001": w_H, "8016": w_O}
        result = fractions_to_iso_densities(comp, "mass_density", rho)

        expected_H = rho * AVOGADRO * w_H / (A_H * CM2_TO_BARN)
        expected_O = rho * AVOGADRO * w_O / (A_O * CM2_TO_BARN)

        assert result["1001"] == pytest.approx(expected_H, rel=1e-10)
        assert result["8016"] == pytest.approx(expected_O, rel=1e-10)

    def test_unknown_density_type_raises(self):
        with pytest.raises(ValueError, match="Unknown density_type"):
            fractions_to_iso_densities({"1001": 1.0}, "unknown", 1.0)

    def test_custom_tolerance(self):
        """With a tight tolerance, 0.999 should be accepted; 0.9 should not."""
        comp = {"1001": 0.999, "8016": 0.001}
        # default tolerance 1e-3 => passes
        fractions_to_iso_densities(comp, "atomic_density", 1.0)

        comp_bad = {"1001": 0.9}
        with pytest.raises(ValueError, match="fractions sum"):
            fractions_to_iso_densities(comp_bad, "atomic_density", 1.0)


# ---------------------------------------------------------------------------
#  parse_all_compositions_from_yaml — integration with GE-14 YAML
# ---------------------------------------------------------------------------
class TestYAMLParsing:
    YAML_PATH = GE14_COMPOSITIONS_YAML

    def test_parse_direct_isotopic(self):
        """UOX fuels have 'isotopic_composition' => parse unchanged."""
        comps = parse_all_compositions_from_yaml(self.YAML_PATH)
        names = [c.material_name for c in comps]
        assert "UOX16" in names
        uox16 = next(c for c in comps if c.material_name == "UOX16")
        # U-235 number density direct from YAML
        assert uox16.isotopic_composition["92235"] == pytest.approx(
            3.758178898925516e-4, rel=1e-10
        )

    def test_parse_atomic_density_moderator(self):
        """MODERATOR uses atomic_density + composition."""
        comps = parse_all_compositions_from_yaml(self.YAML_PATH)
        mod = next(c for c in comps if c.material_name == "MODERATOR")

        N_tot = 7.38974602e-2
        assert mod.isotopic_composition["1001"] == pytest.approx(
            2.0 / 3.0 * N_tot, rel=1e-10
        )
        assert mod.isotopic_composition["8016"] == pytest.approx(
            1.0 / 3.0 * N_tot, rel=1e-10
        )

    def test_depletable_flag(self):
        comps = parse_all_compositions_from_yaml(self.YAML_PATH)
        uox40 = next(c for c in comps if c.material_name == "UOX40")
        mod = next(c for c in comps if c.material_name == "MODERATOR")
        assert uox40.depletable is True
        assert mod.depletable is False

    def test_therm_flag_moderator(self):
        """MODERATOR has ``therm: true`` in the YAML."""
        comps = parse_all_compositions_from_yaml(self.YAML_PATH)
        mod = next(c for c in comps if c.material_name == "MODERATOR")
        assert mod.therm is True
        assert len(mod.therm_data) >= 1
        h1_entry = next(e for e in mod.therm_data if e["isotope"] == "H1")
        assert h1_entry["dragon_alias"] == "H1_H2O"
        assert h1_entry["serpent2_therm_name"] == "lwtr"
        assert h1_entry["serpent2_zaid"] == "1001"

    def test_therm_flag_coolant(self):
        """COOLANT has ``therm: true`` in the YAML."""
        comps = parse_all_compositions_from_yaml(self.YAML_PATH)
        cool = next(c for c in comps if c.material_name == "COOLANT")
        assert cool.therm is True
        assert len(cool.therm_data) >= 1

    def test_therm_flag_fuel_is_false(self):
        """Fuel compositions should NOT have therm set."""
        comps = parse_all_compositions_from_yaml(self.YAML_PATH)
        uox16 = next(c for c in comps if c.material_name == "UOX16")
        assert uox16.therm is False
        assert uox16.therm_data == []


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
