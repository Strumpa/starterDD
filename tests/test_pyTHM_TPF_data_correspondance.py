"""
    Test for validation of correspondance betweenpyTHM (and THM:) and TPF data
    for a given set of input parameters.

    Author.s: A. Mitchell (alec.mitchell@polymtl.ca)
    Created on: 2026-08-19 by A. Mitchell
    Last Modified: 2026-08-19 by A. Mitchell
"""

import os
import re
import pytest

from starterDD.DDModel.DonjonModel import CoreModel
from starterDD.GeometryAnalysis.cartesian_geometry_analysis import CartesianGeometricAnalyser
from starterDD.DDModel.DragonModel import CartesianAssemblyModel

from conftest import (
    GE14_CORE_YAML,
    GE14_SINGLE_ASSEMBLY_CORE_YAML,
    GE14_DOM_GEOMETRY_YAML,
    GE14_VAN_GEOMETRY_YAML
)

# ---------------------------------------------------------------------------
# Test configuration constants
# ---------------------------------------------------------------------------
PATH_TO_YAML_CORE_GEOMETRY = GE14_CORE_YAML
PATH_TO_YAML_SINGLE_ASSEMBLY_CORE_GEOMETRY = GE14_SINGLE_ASSEMBLY_CORE_YAML
PATH_TO_YAML_DOM_ASSEMBLY_GEOMETRY = GE14_DOM_GEOMETRY_YAML
PATH_TO_YAML_VAN_ASSEMBLY_GEOMETRY = GE14_VAN_GEOMETRY_YAML

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def ge14_core():
    """
    Create an object containing the parsed geometry data for the core of the GE14 assembly.
    """
    ge14_core = CoreModel(name="GE14_core", path_to_yaml_configs=PATH_TO_YAML_CORE_GEOMETRY, 
                     core_description_yaml=PATH_TO_YAML_SINGLE_ASSEMBLY_CORE_GEOMETRY)
    return ge14_core

@pytest.fixture(scope="module")
def ge14_core_with_assembly_models(ge14_core):
    """
    Create the assembly models in the core object
    """
    ge14_core.createAssemblyModels()
    return ge14_core

@pytest.fixture(scope="module")
def cartesian_analyser(ge14_core_with_assembly_models):
    """
    Create an object for analyzing the Cartesian geometry of the GE14 core.
    """
    cartesian_analyser = CartesianGeometricAnalyser(core_model=ge14_core_with_assembly_models, core_i=1, core_j=1)
    return cartesian_analyser

@pytest.fixture(scope="module")
def pyTHM_data(cartesian_analyser):
    """
    Create an object containing the pyTHM data with nz=40
    """
    pyTHM_data = cartesian_analyser.run_THM_analysis(40, include_water_rods=False)
    return pyTHM_data

@pytest.fixture(scope="module")
def ge14_cartesian_assembly_models():
    """
    Create the GE14 cartesian assembly model from the DragonModel yaml 
    parser.
    """
    ge14_cartesian_assembly_models = {}
    ge14_cartesian_assembly_models["DOM"] = CartesianAssemblyModel(name="ge14",
                                    tdt_file=None, geometry_description_yaml=PATH_TO_YAML_DOM_ASSEMBLY_GEOMETRY)
    ge14_cartesian_assembly_models["VAN"] = CartesianAssemblyModel(name="ge14",
                                    tdt_file=None, geometry_description_yaml=PATH_TO_YAML_VAN_ASSEMBLY_GEOMETRY)
    return ge14_cartesian_assembly_models

@pytest.fixture(scope="module")
def ge14_cartesian_assembly_models_parsed(ge14_cartesian_assembly_models):
    """Parse the GE14 DOM and VAN assembly models."""
    ge14_cartesian_assembly_models["DOM"].parse_geometry_description(PATH_TO_YAML_DOM_ASSEMBLY_GEOMETRY)
    ge14_cartesian_assembly_models["VAN"].parse_geometry_description(PATH_TO_YAML_VAN_ASSEMBLY_GEOMETRY)
    ge14_cartesian_assembly_models_parsed = {
        "DOM": ge14_cartesian_assembly_models["DOM"],
        "VAN": ge14_cartesian_assembly_models["VAN"]
    }
    return ge14_cartesian_assembly_models_parsed

@pytest.fixture(scope="module")
def tpf_data(cartesian_analyser):
    """
    Create an object containing the TPF data with nz=40
    """
    tpf_data = cartesian_analyser.run_TPF_analysis(40)
    return tpf_data

# ----------------------------------------------------------------------------
# Tests: pyTHM geometry processing
# ----------------------------------------------------------------------------
class TestpyTHMGeometry:
    """Tests for pyTHM geometry processing"""
    def test_pyTHM_geometric_data_exists(self, pyTHM_data):
        """Verify that pyTHM geometric data exists"""
        assert pyTHM_data is not None

    def test_pyTHM_fuel_data(self, pyTHM_data, ge14_cartesian_assembly_models_parsed):
        """Verify that pyTHM fuel data corresponds to
        the expected geometry."""
        ge14_dom_assembly_parsed = ge14_cartesian_assembly_models_parsed["DOM"]
        expected_pin_pitch = ge14_dom_assembly_parsed.pin_geometry_dict.get("pin_pitch", 0) * 1E-2
        assert pyTHM_data["fuel_data"]["pin_pitch"] == expected_pin_pitch
        expected_gap_radius = ge14_dom_assembly_parsed.pin_geometry_dict.get("gap_radius", None) * 1E-2
        assert pyTHM_data["fuel_data"]["gap_radius"] == expected_gap_radius

# ----------------------------------------------------------------------------
# Tests: pyTHM - TPF comparison
# ----------------------------------------------------------------------------
class TestpyTHMTPF:
    """Tests for comparing pyTHM and TPF data"""
    def test_tpf_data_exists(self, tpf_data):
        """Verify that TPF data exists"""
        assert tpf_data is not None
    def test_pyTHM_tpf_comparison(self, tpf_data, pyTHM_data):
        """Verify that pyTHM and TPF data correspond to each other"""
        assert pyTHM_data["active_flow_data"]["pitch"] == tpf_data["active_flow_data"]["pitch"]
        assert pyTHM_data["fuel_data"]["pin_pitch"] == tpf_data["fuel_data"]["pin_pitch"]
        assert pyTHM_data["fuel_data"]["gap_radius"] == tpf_data["fuel_data"]["gap_radius"]
        expected_height = 347.1291 # Taken directly from the yaml instead of parsing
        assert tpf_data["active_flow_data"]["core_height"] == expected_height * 1E-2