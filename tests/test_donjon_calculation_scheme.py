import pytest 
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from starterDD.DDModel.DonjonModel import CoreModel
from conftest import GE14_CORE_YAML
from starterDD.InterfaceToDD.donjon_module_calls import INIT, NEUTRONICS
from starterDD.DDModel.DonjonCalculationScheme import DonjonCalculationScheme, ModelInitialisation, NeutronicsSolve


@pytest.fixture
def one_channel_core_model():
    GE14_core_description_yaml = "GEOM_single_assembly_CORE.yaml"
    one_channel_core_model = CoreModel(name="GE14_single_assembly_core", path_to_yaml_configs=GE14_CORE_YAML, core_description_yaml=GE14_core_description_yaml)
    return one_channel_core_model

@pytest.fixture
def DOM_channel_initialiser():
    GE14_core_description_yaml = "GEOM_single_DOM_assembly_CORE.yaml"
    core_model = CoreModel(name="GE14_single_DOM_assembly_core", path_to_yaml_configs=GE14_CORE_YAML, core_description_yaml=GE14_core_description_yaml)
    DOM_channel_initialiser = ModelInitialisation(core_model=core_model, 
                                        radial_homogenization_strategy="by_assembly",
                                        boundary_conditions_dict={"x": "reflective",
                                                                  "y": "reflective",
                                                                  "z": "void"},
                                        number_of_axial_materials_per_slice=[40],
                                        number_of_energy_groups=2,
                                        interpolation_variables=[("TFuel", "local", "fuel_temperature"), 
                                                                 ("TCool", "local", "coolant_temperature"),
                                                                 ("DCool", "local", "coolant_density")])
    
    DOM_channel_initialiser.resolve_mesh()

    DOM_channel_initialiser.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[0.0, 347.1291], compo_name="CPODOM", dir_name="EDI_2G")

    DOM_channel_initialiser.set_reactor_power(3.6e6)
    DOM_channel_initialiser.set_initial_axial_power_form("uniform")
    DOM_channel_initialiser.set_initial_parameters({"TFuel": 900, "TCool": 600.0, "DCool": 0.73669})
    DOM_channel_initialiser.set_fuel_mass(total_fuel_mass = 0.5)
    
    return DOM_channel_initialiser

@pytest.fixture
def model_initialiser():
    GE14_core_description_yaml = "GEOM_single_assembly_CORE.yaml"
    one_channel_core_model = CoreModel(name="GE14_single_assembly_core", path_to_yaml_configs=GE14_CORE_YAML, core_description_yaml=GE14_core_description_yaml)
    model_initialiser = ModelInitialisation(core_model=one_channel_core_model, 
                                        radial_homogenization_strategy="by_pin",
                                        boundary_conditions_dict={"x": "reflective",
                                                                  "y": "reflective",
                                                                  "z": "void"},
                                        number_of_axial_materials_per_slice=[10, 5],
                                        number_of_energy_groups=2,
                                        interpolation_variables=[("local", "TFuel"), 
                                                                 ("local", "TCool"),
                                                                 ("local", "DCool")])
    
    model_initialiser.resolve_mesh()

    model_initialiser.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[0.0, 222.0595], compo_name="CPODOM", dir_name="EDI_2G")
    model_initialiser.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[222.0595, 347.1291], compo_name="CPOVAN", dir_name="EDI_2G")
    
    return model_initialiser


@pytest.fixture
def minicore_initialiser():
    GE14_core_description_yaml = "GEOM_mini_CORE_reflector.yaml"
    mini_core_model = CoreModel(name="GE14_4x4_core_reflectors", path_to_yaml_configs=GE14_CORE_YAML, core_description_yaml=GE14_core_description_yaml)
    minicore_initialiser = ModelInitialisation(core_model=mini_core_model, 
                                    radial_homogenization_strategy="by_assembly",
                                    boundary_conditions_dict={"x": "reflective",
                                                              "y": "reflective",
                                                              "z": "void"},
                                    number_of_axial_materials_per_slice=[10, 5],
                                    number_of_energy_groups=2,
                                    interpolation_variables=[("TFuel", "local", "fuel_temperature"), 
                                                                ("TCool", "local", "coolant_temperature"),
                                                                ("DCool", "local", "coolant_density")])
    minicore_initialiser.resolve_mesh()
    minicore_initialiser.set_slice_to_compodir_correspondance(xpos="ALL", ypos="ALL", z_bounds=[0.0, 222.0595], compo_name="CPODOM", dir_name="EDI_2G")
    minicore_initialiser.set_slice_to_compodir_correspondance(xpos="ALL", ypos="ALL", z_bounds=[222.0595, 347.1291], compo_name="CPOVAN", dir_name="EDI_2G")
    
    return minicore_initialiser


def generate_macro_mesh(xmesh, ymesh, split_diag, tol=1e-12):
    """
    Generate the mesh elements lying below the symmetry axis.

    The symmetry axis is the line connecting:
        (xmin, ymin) -> (xmax, ymax)

    Elements are returned as:

        RECT x1 x2 y1 y2

    or

        TRIA x1 y1 x2 y2 x3 y3

    Only the portion below the symmetry axis is retained.

    Parameters
    ----------
    xmesh : array-like
        Mesh points in x.
    ymesh : array-like
        Mesh points in y.
    tol : float
        Numerical tolerance.

    Returns
    -------
    elements : list[str]
        Mesh elements below the symmetry axis.
    """

    xmin, xmax = xmesh[0], xmesh[-1]
    ymin, ymax = ymesh[0], ymesh[-1]

    # Slope of the symmetry axis
    slope = (ymax - ymin) / (xmax - xmin)

    def f(x, y):
        """
        Signed value relative to the symmetry line.

        f < 0 : below the symmetry axis
        f = 0 : on the symmetry axis
        f > 0 : above the symmetry axis
        """
        y_line = ymin + slope * (x - xmin)
        return y - y_line

    def line_intersection(P, Q):
        """
        Find the intersection between segment P-Q and the
        symmetry line.
        """

        fP = f(*P)
        fQ = f(*Q)

        t = fP / (fP - fQ)

        return (
            P[0] + t * (Q[0] - P[0]),
            P[1] + t * (Q[1] - P[1])
        )

    def clip_below_line(polygon):
        """
        Clip polygon against the half-plane below the
        symmetry line using Sutherland-Hodgman clipping.
        """
        result = []
        for i in range(len(polygon)):
            P = polygon[i]
            Q = polygon[(i + 1) % len(polygon)]
            fP = f(*P)
            fQ = f(*Q)
            inside_P = fP <= tol
            inside_Q = fQ <= tol

            if inside_P:
                result.append(P)

            # Edge crosses the symmetry line
            if inside_P != inside_Q:
                result.append(line_intersection(P, Q))
        return result

    elements = []
    for i in range(len(xmesh) - 1):
        x1 = xmesh[i]
        x2 = xmesh[i + 1]
        for j in range(len(ymesh) - 1):
            y1 = ymesh[j]
            y2 = ymesh[j + 1]
            rectangle = [(x1, y1),(x2, y1),(x2, y2),(x1, y2)]
            values = [f(x, y) for x, y in rectangle]

            # ---------------------------------------------------------
            # Case 1: Entire rectangle is above the symmetry axis
            # ---------------------------------------------------------
            if split_diag:
                if all(value > tol for value in values):
                    continue

            if split_diag:
                # ---------------------------------------------------------
                # Case 2: Entire rectangle is below the symmetry axis
                # ---------------------------------------------------------
                if all(value <= tol for value in values):
                    elements.append(f"RECT {x1} {x2} {y1} {y2}")
                    continue

                # ---------------------------------------------------------
                # Case 3: Symmetry axis cuts the rectangle
                # ---------------------------------------------------------

                polygon = clip_below_line(rectangle)
                # Remove duplicate points
                cleaned = []
                for point in polygon:
                    if not cleaned:
                        cleaned.append(point)
                        continue
                    if (abs(point[0] - cleaned[-1][0]) > tol or abs(point[1] - cleaned[-1][1]) > tol):
                        cleaned.append(point)

                # Remove closing duplicate
                if len(cleaned) > 1:
                    if (abs(cleaned[0][0] - cleaned[-1][0]) <= tol and abs(cleaned[0][1] - cleaned[-1][1]) <= tol):
                        cleaned.pop()

                # ---------------------------------------------------------
                # The retained portion is a triangle
                # ---------------------------------------------------------
                if len(cleaned) == 3:
                    p1, p2, p3 = cleaned
                    elements.append(
                        f"TRIA {p1[0]} {p1[1]} {p2[0]} {p2[1]} {p3[0]} {p3[1]}")

                # ---------------------------------------------------------
                # The retained portion is a quadrilateral.
                # Split it into two triangles.
                # ---------------------------------------------------------
                elif len(cleaned) == 4:
                    p1, p2, p3, p4 = cleaned
                    elements.append(
                        f"TRIA {p1[0]} {p1[1]} {p2[0]} {p2[1]} {p3[0]} {p3[1]}")
                    elements.append(
                        f"TRIA {p1[0]} {p1[1]} {p3[0]} {p3[1]} {p4[0]} {p4[1]}")
            else:
                elements.append(f"RECT {x1} {x2} {y1} {y2}")

    print(len(elements))
    cleaned_elements = "\n".join(elements)
    return cleaned_elements, elements


def test_model_initialisation(one_channel_core_model):
    # Initialisation step creation
    init_step = ModelInitialisation(core_model=one_channel_core_model, 
                                    radial_homogenization_strategy="by_assembly",
                                    boundary_conditions_dict={"x": "reflective",
                                                              "y": "reflective",
                                                              "z": "void"},
                                    number_of_axial_materials_per_slice=[10, 5],
                                    number_of_energy_groups=2,
                                    interpolation_variables=[("TFuel", "local", "fuel_temperature"), 
                                                            ("TCool", "local", "coolant_temperature"),
                                                            ("DCool", "local", "coolant_density")])

    init_step.resolve_mesh()

    assert init_step.meshx == [0.0, 15.24]
    assert init_step.meshy == [0.0, 15.24]
    assert len(init_step.meshz) == 16

    init_step = ModelInitialisation(core_model=one_channel_core_model, 
                                    radial_homogenization_strategy="by_pin",
                                    boundary_conditions_dict={"x": "reflective",
                                                              "y": "reflective",
                                                              "z": "void"},
                                    number_of_axial_materials_per_slice=[20, 30],
                                    number_of_energy_groups=2,
                                    interpolation_variables=[("TFuel", "local", "fuel_temperature"), 
                                                            ("TCool", "local", "coolant_temperature"),
                                                            ("DCool", "local", "coolant_density")])

    init_step.resolve_mesh()
    ref_mesh = [0.0, 1.12, 2.42, 3.72, 5.02, 6.32, 7.62, 8.92, 10.22, 11.52, 12.82, 14.12, 15.24]
    for meshpt_x, refmeshpt in zip(init_step.meshx, ref_mesh):
        assert meshpt_x == pytest.approx(refmeshpt, abs=1e-5) 
    for meshpt_y, refmeshpt in zip(init_step.meshy, ref_mesh):
        assert meshpt_y == pytest.approx(refmeshpt, abs=1e-5) 
    assert len(init_step.meshz) == 51

    init_step.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[0.0, 222.0595], compo_name="CPODOM", dir_name="EDI_2G")
    init_step.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[222.0595, 347.1291], compo_name="CPOVAN", dir_name="EDI_2G")

    elements, raw_elements = generate_macro_mesh(ref_mesh, ref_mesh, split_diag=True)
    init_step.set_reactor_power(870e6*16/240)


def test_core_model_analysis(minicore_initialiser):
    assert len(minicore_initialiser.meshz) == 16
    assert len(minicore_initialiser.meshy) == 7
    assert len(minicore_initialiser.meshx) == 7
    ref_meshx = [0.0, 15.24, 30.48, 45.72, 60.96, 76.2, 91.44]
    for meshpt, refpt in zip(minicore_initialiser.meshx, ref_meshx):
        assert meshpt == pytest.approx(refpt, abs=1e-5)

    assert len(minicore_initialiser.reflector_pairs) == 6*6 - 4*4
    assert len(minicore_initialiser.fuel_channel_pairs) == 4*4


def test_neutronics_solve(minicore_initialiser):
    neutonics_solution = NeutronicsSolve(minicore_initialiser, "diffusion", "linear", ["MCFD", 1])

    assert neutonics_solution.interpolation_type == "linear"
    assert neutonics_solution.operator == "diffusion"

    nbFuelChannels = minicore_initialiser.number_fuel_channels
    nz = len(minicore_initialiser.meshz) - 1
    assert "TFuel" in neutonics_solution.model.interpolation_variables.keys()
    assert neutonics_solution.model.interpolation_variables["TFuel"][0] == "LOCAL"
    neutonics_solution.set_interpolation_parameter("TFuel", [1200]*nbFuelChannels*nz)
    assert len(neutonics_solution.parameter_values_dict["TFuel"]) == nbFuelChannels*nz
    neutonics_solution.resolve_cpo_dir_to_mix()
    assert len(neutonics_solution.compo_dir_to_mix) == 2
    assert len(neutonics_solution.compo_dir_to_mix[('CPODOM', 'EDI_2G')]) == 160
    assert len(neutonics_solution.compo_dir_to_mix[('CPOVAN', 'EDI_2G')]) == 240 - 160

def test_single_DOM_channel(DOM_channel_initialiser):


    neutronics_step = NeutronicsSolve(DOM_channel_initialiser, "diffusion", "linear", ["MCFD", 1])
    nz = len(DOM_channel_initialiser.meshz) - 1
    neutronics_step.set_interpolation_parameter("TFuel", [1200.0]*nz)
    neutronics_step.set_interpolation_parameter("TCool", [559.0]*nz)
    neutronics_step.set_interpolation_parameter("DCool", [0.600]*nz)
    neutronics_step.resolve_cpo_dir_to_mix()

    scheme = DonjonCalculationScheme(name="test_scheme")
    scheme.add_initialisation_step(DOM_channel_initialiser)
    scheme.add_neutronics_step(neutronics_step)

    init_proc = INIT(scheme)
    assert init_proc.scheme.initialisation_step.nfuel == 40

    neutronics_proc = NEUTRONICS(scheme)
    neutronics_proc.build_RESINI_call()
    neutronics_proc.build_NCR_call()

    init_proc.write_to_c2m("tests/outputs", "IniDonjon_minicore")
