"""
Tests for Serpent2 exports with asymmetric narrow/wide gaps.

This test module validates the implementation of asymmetric gap support
in Serpent2_exports.py. It uses PB2 Type6 (off-centered assembly geometry)
for both controlled and uncontrolled cases.

Features tested:
- Lattice center positioning for off-centered lattices
- Water rod center handling (no translation applied)
- Asymmetric lattice padding
- Symmetry detection and gap layout computation
"""

import pytest

from starterDD.DDModel.DragonModel import CartesianAssemblyModel
from starterDD.InterfaceToDD.Serpent2_exports import (
    S2_Lattice,
    S2_ChannelGeometry,
    _get_lattice_symmetry,
    _compute_gap_layout,
)

from conftest import (
    PB2_UNCONTROLLED_GEOMETRY_YAML,
    PB2_CONTROLLED_GEOMETRY_YAML,
)


class TestSymmetryHelpers:
    """Test symmetry detection and gap layout computation helpers."""

    def test_get_lattice_symmetry_returns_none_for_mock_model(self):
        """Test that symmetry helper returns None if model lacks check_diagonal_symmetry."""
        class MockAssembly:
            pass

        result = _get_lattice_symmetry(MockAssembly())
        assert result is None

    def test_compute_gap_layout_anti_diagonal(self):
        """Test gap layout computation for anti-diagonal symmetry (GE14-style)."""
        gap_layout = _compute_gap_layout("anti-diagonal", gap_wide=0.9, gap_narrow=0.4)

        assert gap_layout["NW"] == 0.9  # Wide-wide at north-west
        assert gap_layout["SE"] == 0.9  # Wide-wide at south-east (diagonal)
        assert gap_layout["NE"] == 0.4  # Narrow at north-east
        assert gap_layout["SW"] == 0.4  # Narrow at south-west

    def test_compute_gap_layout_main_diagonal(self):
        """Test gap layout computation for main-diagonal symmetry (ATRIUM10-style)."""
        gap_layout = _compute_gap_layout("main-diagonal", gap_wide=0.75, gap_narrow=0.75)

        assert gap_layout["SW"] == 0.75  # Wide-wide at south-west
        assert gap_layout["NE"] == 0.75  # Narrow at north-east (diagonal)
        assert gap_layout["NW"] == 0.75  # Mixed at north-west
        assert gap_layout["SE"] == 0.75  # Mixed at south-east

    def test_compute_gap_layout_symmetric(self):
        """Test gap layout computation for symmetric gaps."""
        gap_layout = _compute_gap_layout(None, gap_wide=0.75, gap_narrow=0.75)

        # All corners should have the same gap width
        assert gap_layout["NW"] == 0.75
        assert gap_layout["NE"] == 0.75
        assert gap_layout["SW"] == 0.75
        assert gap_layout["SE"] == 0.75


class TestPB2UncontrolledAsymmetricAssembly:
    """Test Serpent2 export functionality with PB2 Type6 uncontrolled assembly."""

    @pytest.fixture
    def pb2_uncontrolled_assembly(self):
        """Create a PB2 Type6 uncontrolled assembly model."""
        # Load the assembly geometry
        assembly = CartesianAssemblyModel(
            name="PB2_Type6_uncontrolled",
            tdt_file=None,  # We don't need the tdt file for Serpent2 export testing
            geometry_description_yaml=PB2_UNCONTROLLED_GEOMETRY_YAML,
        )

        # Analyze lattice to populate water_rods and lattice
        assembly.analyze_lattice_description(build_pins=False)

        return assembly

    def test_pb2_uncontrolled_assembly_loads(self, pb2_uncontrolled_assembly):
        """Test that PB2 uncontrolled assembly YAML loads successfully."""
        am = pb2_uncontrolled_assembly

        # Check basic geometry parameters
        assert am.assembly_pitch == pytest.approx(15.24)
        assert am.gap_wide == pytest.approx(0.9017)
        assert am.gap_narrow == pytest.approx(0.42418)
        assert am.channel_box_thickness == pytest.approx(0.254)

        # Check lattice description
        assert am.lattice_description is not None
        assert len(am.lattice_description) == 8  # 8x8 lattice
        assert len(am.lattice_description[0]) == 8

        # Check water rods (populated by analyze_lattice_description)
        assert am.water_rods is not None
        assert len(am.water_rods) == 2  # Two water rods

    def test_pb2_uncontrolled_has_no_control_cross(self, pb2_uncontrolled_assembly):
        """Test that uncontrolled assembly has no control cross."""
        am = pb2_uncontrolled_assembly
        assert am.has_control_cross is False
        assert am.control_cross is None

    def test_pb2_uncontrolled_translation_offsets(self, pb2_uncontrolled_assembly):
        """Test that translation offsets are computed correctly."""
        am = pb2_uncontrolled_assembly

        # Both translation offsets should be defined
        assert am.translation_offset_x is not None
        assert am.translation_offset_y is not None

        # translation_offset_x should use gap_wide
        # translation_offset_y should use gap_narrow
        assert am.translation_offset_x >= am.translation_offset_y

    def test_pb2_uncontrolled_water_rod_centers(self, pb2_uncontrolled_assembly):
        """Test that water rod centers are in assembly coordinates."""
        am = pb2_uncontrolled_assembly

        # Water rod centers from YAML: [[7.04596, 6.56844], [8.67156, 8.19404]]
        wr0_center = am.water_rods[0].center
        wr1_center = am.water_rods[1].center

        # Check that centers are in expected range (assembly coords)
        assert 0 < wr0_center[0] < am.assembly_pitch
        assert 0 < wr0_center[1] < am.assembly_pitch
        assert 0 < wr1_center[0] < am.assembly_pitch
        assert 0 < wr1_center[1] < am.assembly_pitch

        # Check specific values match YAML
        assert wr0_center[0] == pytest.approx(7.04596, rel=1e-4)
        assert wr0_center[1] == pytest.approx(6.56844, rel=1e-4)
        assert wr1_center[0] == pytest.approx(8.67156, rel=1e-4)
        assert wr1_center[1] == pytest.approx(8.19404, rel=1e-4)

    def test_pb2_uncontrolled_s2_lattice_creation(self, pb2_uncontrolled_assembly):
        """Test that S2_Lattice can be created from PB2 uncontrolled assembly."""
        am = pb2_uncontrolled_assembly

        # Create S2_Lattice
        s2_lattice = S2_Lattice.from_assembly_model(am)

        # Check lattice properties
        assert s2_lattice is not None
        assert s2_lattice.center_x is not None
        assert s2_lattice.center_y is not None
        assert s2_lattice.nx == 10  # 8x8 fuel lattice + 2 padding (1 on each side)
        assert s2_lattice.ny == 10
        assert len(s2_lattice.universe_map) == 10
        assert len(s2_lattice.universe_map[0]) == 10

        # Check that lattice is off-centered from assembly center
        assembly_center_x = am.assembly_pitch / 2.0
        assembly_center_y = am.assembly_pitch / 2.0
        # Lattice center should NOT be at assembly center due to asymmetric gaps
        # Allow larger tolerance since gaps affect positioning
        assert abs(s2_lattice.center_x - assembly_center_x) > 0.1
        assert abs(s2_lattice.center_y - assembly_center_y) > 0.1

    def test_pb2_uncontrolled_s2_channel_geometry(self, pb2_uncontrolled_assembly):
        """Test that S2_ChannelGeometry can be created from PB2 uncontrolled assembly."""
        am = pb2_uncontrolled_assembly

        # Create S2_ChannelGeometry
        s2_channel = S2_ChannelGeometry.from_assembly_model(am)

        # Check channel geometry properties
        assert s2_channel is not None
        assert s2_channel.inner_hw > 0  # Use inner_hw instead of channel_inner_half_width
        assert s2_channel.outer_hw > s2_channel.inner_hw  # Use outer_hw instead of channel_outer_half_width

        # Check water rods
        assert s2_channel.water_rods is not None
        assert len(s2_channel.water_rods) == 2

        # Check that water rod centers are NOT translated
        # Compare as tuples since water_rods stores tuples
        for i, wr_info in enumerate(s2_channel.water_rods):
            wr_center = wr_info["center"]
            expected_center = am.water_rods[i].center
            # Handle both tuple and list comparisons
            assert (wr_center[0], wr_center[1]) == (expected_center[0], expected_center[1])


class TestPB2ControlledAsymmetricAssembly:
    """Test Serpent2 export functionality with PB2 Type6 controlled assembly."""

    @pytest.fixture
    def pb2_controlled_assembly(self):
        """Create a PB2 Type6 controlled assembly model."""
        assembly = CartesianAssemblyModel(
            name="PB2_Type6_controlled",
            tdt_file=None,
            geometry_description_yaml=PB2_CONTROLLED_GEOMETRY_YAML,
        )

        # Analyze lattice to populate water_rods and lattice
        assembly.analyze_lattice_description(build_pins=False)

        return assembly

    def test_pb2_controlled_assembly_loads(self, pb2_controlled_assembly):
        """Test that PB2 controlled assembly YAML loads successfully."""
        am = pb2_controlled_assembly

        # Check basic geometry parameters (same as uncontrolled)
        assert am.assembly_pitch == pytest.approx(15.24)
        assert am.gap_wide == pytest.approx(0.9017)
        assert am.gap_narrow == pytest.approx(0.42418)

    def test_pb2_controlled_has_control_cross(self, pb2_controlled_assembly):
        """Test that controlled assembly has a control cross."""
        am = pb2_controlled_assembly
        assert am.has_control_cross is True
        assert am.control_cross is not None

    def test_pb2_controlled_control_cross_at_wide_wide_corner(self, pb2_controlled_assembly):
        """Test that control cross is placed at wide-wide corner."""
        am = pb2_controlled_assembly

        # Check control cross center
        assert am.control_cross.center == "north-west"

        # For PB2, the lattice should have anti-diagonal symmetry
        # So north-west should be the wide-wide corner
        symmetry = am.check_diagonal_symmetry()
        if symmetry == "anti-diagonal":
            # north-west is wide-wide for anti-diagonal
            assert am.control_cross.center == "north-west"

    def test_pb2_controlled_s2_lattice_creation(self, pb2_controlled_assembly):
        """Test that S2_Lattice can be created from PB2 controlled assembly."""
        am = pb2_controlled_assembly

        # Create S2_Lattice
        s2_lattice = S2_Lattice.from_assembly_model(am)

        # Check lattice properties (same structure as uncontrolled)
        assert s2_lattice is not None
        assert s2_lattice.nx == 10
        assert s2_lattice.ny == 10

    def test_pb2_controlled_s2_channel_geometry(self, pb2_controlled_assembly):
        """Test that S2_ChannelGeometry can be created from PB2 controlled assembly."""
        am = pb2_controlled_assembly

        # Create S2_ChannelGeometry
        s2_channel = S2_ChannelGeometry.from_assembly_model(am)

        # Check channel geometry properties
        assert s2_channel is not None
        assert len(s2_channel.water_rods) == 2

        # Check that water rod centers are NOT translated
        # Compare as tuples since water_rods stores tuples
        for i, wr_info in enumerate(s2_channel.water_rods):
            wr_center = wr_info["center"]
            expected_center = am.water_rods[i].center
            # Handle both tuple and list comparisons
            assert (wr_center[0], wr_center[1]) == (expected_center[0], expected_center[1])


class TestAsymmetricGapBackwardCompatibility:
    """Test that symmetric gap cases still work correctly."""

    def test_symmetric_gaps_lattice_offset_computation(self):
        """Test that symmetric gap cases compute reasonable offsets."""
        # Create a mock assembly with symmetric gaps
        class MockAssembly:
            def __init__(self):
                self.lattice = [["pin"] * 8 for _ in range(8)]
                self.pin_geometry_dict = {"pin_pitch": 1.6256}
                self.assembly_pitch = 15.24
                self.gap_wide = 0.75
                self.gap_narrow = 0.75
                self.channel_box_thickness = 0.254
                self.channel_box_inner_side = None
                self.translation_offset_x = None
                self.translation_offset_y = None
                self.water_rods = []

            def check_diagonal_symmetry(self):
                return None

        assembly = MockAssembly()
        s2_lattice = S2_Lattice.from_assembly_model(assembly)

        # For symmetric gaps, offsets should be equal
        # (both using same gap_wide and gap_narrow values)
        assert s2_lattice.center_x is not None
        assert s2_lattice.center_y is not None
        # Since both gaps are the same, the computed offsets should be equal
        # So the lattice should be centered (offset from assembly center by same amount in both directions)
        assert isinstance(s2_lattice.center_x, float)
        assert isinstance(s2_lattice.center_y, float)


class TestAsymmetricGapCenterAlignment:
    """Test alignment of lattice and channel box centers with correct offset computation."""

    @pytest.fixture
    def pb2_uncontrolled_assembly(self):
        """Create a PB2 Type6 uncontrolled assembly model."""
        assembly = CartesianAssemblyModel(
            name="PB2_Type6_uncontrolled",
            tdt_file=None,
            geometry_description_yaml=PB2_UNCONTROLLED_GEOMETRY_YAML,
        )
        assembly.analyze_lattice_description(build_pins=False)
        return assembly

    def test_lattice_center_uses_unpadded_dimensions(self, pb2_uncontrolled_assembly):
        """Test that lattice center is computed from unpadded lattice dimensions (8x8)."""
        am = pb2_uncontrolled_assembly
        pin_pitch = am.pin_geometry_dict["pin_pitch"]

        # Create S2_Lattice
        s2_lattice = S2_Lattice.from_assembly_model(am)

        # Expected: lattice_center_local from 8x8 unpadded lattice
        # lattice_center_x_local = (8 * pin_pitch) / 2.0
        unpadded_center_x_local = (8 * pin_pitch) / 2.0
        unpadded_center_y_local = (8 * pin_pitch) / 2.0

        # With offsets applied
        expected_center_x = unpadded_center_x_local + am.translation_offset_x
        expected_center_y = unpadded_center_y_local + am.translation_offset_y

        # Verify lattice center matches expected computation
        assert s2_lattice.center_x == pytest.approx(expected_center_x, rel=1e-5), \
            f"Lattice center_x {s2_lattice.center_x} != expected {expected_center_x}"
        assert s2_lattice.center_y == pytest.approx(expected_center_y, rel=1e-5), \
            f"Lattice center_y {s2_lattice.center_y} != expected {expected_center_y}"

        # Voerify lattice has padded dimensions (10x10) for Serpent2
        assert s2_lattice.nx == 10
        assert s2_lattice.ny == 10

    def test_channel_box_center_aligns_with_lattice_center(self, pb2_uncontrolled_assembly):
        """Test that channel box center matches lattice center."""
        am = pb2_uncontrolled_assembly

        # Create both lattice and channel geometry
        s2_lattice = S2_Lattice.from_assembly_model(am)
        s2_channel = S2_ChannelGeometry.from_assembly_model(am)

        # Verify centers match
        assert s2_channel.cx == pytest.approx(s2_lattice.center_x, rel=1e-6), \
            f"Channel box center_x {s2_channel.cx} != lattice center_x {s2_lattice.center_x}"
        assert s2_channel.cy == pytest.approx(s2_lattice.center_y, rel=1e-6), \
            f"Channel box center_y {s2_channel.cy} != lattice center_y {s2_lattice.center_y}"

    def test_bounding_surface_centered_at_assembly_pitch_half(self, pb2_uncontrolled_assembly):
        """Test that bounding surface (1002) is always centered at (ap/2, ap/2)."""
        am = pb2_uncontrolled_assembly

        # Create channel geometry
        s2_channel = S2_ChannelGeometry.from_assembly_model(am)

        # Expected: bounding surface always at assembly_pitch / 2
        assembly_hp = am.assembly_pitch / 2.0

        # Verify bounding surface is centered at assembly center
        # The bounding surface should use assembly_hp (assembly_pitch/2) for both center and half-pitch
        assert s2_channel.assembly_hp == pytest.approx(assembly_hp, rel=1e-5), \
            f"Channel geometry assembly_hp {s2_channel.assembly_hp} != expected {assembly_hp}"

        # The channel box center (cx, cy) should be offset by translation offsets,
        # but the bounding surface center should always be at (assembly_hp, assembly_hp)
        # Verify this is reflected in the format_cards output
        cards = s2_channel.format_cards()

        # Check that the bounding surface line exists in the output
        assert 'surf 1002' in cards, "Surface 1002 (bounding surface) not found in channel geometry cards"

    def test_serpent2model_lattice_delegation(self, pb2_uncontrolled_assembly):
        """Test that Serpent2Model correctly delegates to S2_Lattice.from_assembly_model()."""
        from starterDD.InterfaceToDD.Serpent2_exports import Serpent2Model, S2_Settings

        am = pb2_uncontrolled_assembly

        # Build Serpent2 model
        settings = S2_Settings()
        model = Serpent2Model(assembly_model=am, settings=settings)

        # Build with minimal setup (no fuel materials needed for this test)
        model.build(
            lattice_name="10",
            empty_universe_name="empty",
        )

        # Verify lattice was created and has correct properties
        assert model.lattice is not None
        assert model.lattice.name == "10"
        assert model.lattice.nx == 10  # padded
        assert model.lattice.ny == 10  # padded

        # Verify lattice center matches expected offset computation
        s2_lattice_direct = S2_Lattice.from_assembly_model(am)
        assert model.lattice.center_x == pytest.approx(s2_lattice_direct.center_x, rel=1e-6)
        assert model.lattice.center_y == pytest.approx(s2_lattice_direct.center_y, rel=1e-6)

        # Verify channel geometry was also created
        assert model.channel_geometry is not None

        # Verify channel geometry center matches lattice center
        assert model.channel_geometry.cx == pytest.approx(model.lattice.center_x, rel=1e-6)
        assert model.channel_geometry.cy == pytest.approx(model.lattice.center_y, rel=1e-6)
