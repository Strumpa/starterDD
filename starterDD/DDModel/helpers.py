# helpers to build Dragon Models.
# R.Guasch
# Date : 09/02/2026 [created]

import yaml

def associate_material_to_rod_ID(path_to_materials_yaml, path_to_geometry_yaml):
    """
    Associate materials to rod IDs when using a material composition definition with the rod_id entry.
    This is intended to facilitate Models based on benchmark problems where the material composition definition is based on the rod IDs, which is a common approach for BWR benchmarks.
    """

    file_materials = open(path_to_materials_yaml, 'r')
    materials_data = yaml.safe_load(file_materials)
    file_materials.close()
    file_geometry = open(path_to_geometry_yaml, 'r')
    geometry_data = yaml.safe_load(file_geometry)
    file_geometry.close()  

    rod_id_to_material = {}
    mix_list = materials_data.get('MIX_COMPOSITIONS', [])
    for entry in mix_list:
        name = entry.get('name')
        rod_id = entry.get('rod_id')
        if rod_id is not None:
            rod_id_to_material[rod_id] = name 

    # check that all rod IDs in the geometry description have an associated material
    lattice_desc = geometry_data.get('lattice_description', [])
    for row in lattice_desc:
        for rod_id in row:
            if rod_id not in rod_id_to_material:
                raise ValueError(f"No material associated to rod ID {rod_id} found in material compositions YAML file.")
            
    return rod_id_to_material




