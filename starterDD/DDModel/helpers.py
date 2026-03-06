# helpers to build Dragon Models.
# R.Guasch
# Date : 09/02/2026 [created]

import yaml

def associate_material_to_rod_ID(path_to_materials_yaml, path_to_geometry_yaml):
    """Map lattice rod IDs to material names from YAML definitions.

    Reads the ``MIX_COMPOSITIONS`` section of *path_to_materials_yaml*
    for entries that carry a ``rod_id`` key, and builds a
    ``{rod_id: material_name}`` mapping.  Then validates that every
    rod ID appearing in the ``lattice_description`` of
    *path_to_geometry_yaml* has an associated material.

    Parameters
    ----------
    path_to_materials_yaml : str
        Path to the YAML file containing ``MIX_COMPOSITIONS`` entries
        with ``name`` and ``rod_id`` keys.
    path_to_geometry_yaml : str
        Path to the geometry YAML file containing
        ``lattice_description`` (2-D list of rod ID strings).

    Returns
    -------
    dict
        ``{rod_id: material_name}``

    Raises
    ------
    ValueError
        If a rod ID in the lattice has no corresponding material, or
        if either YAML file cannot be read.
    """
    try:
        file_materials = open(path_to_materials_yaml, 'r')
        materials_data = yaml.safe_load(file_materials)
        file_materials.close()
    except Exception as e:
        raise ValueError(f"Error reading materials YAML file: {e}")
    
    try:
        file_geometry = open(path_to_geometry_yaml, 'r')
        geometry_data = yaml.safe_load(file_geometry)
        file_geometry.close()  
    except Exception as e:
        raise ValueError(f"Error reading geometry YAML file: {e}")

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




