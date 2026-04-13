### Collection of classes for basic surface creation in the geometry builder

class Square:
    def __init__(self, name, surface_id, center, side_length, corner_radius_of_curvature=None):
        """
        Define a square surface with : 
        - name : user friendly name for the surface : 
            Rescritions : "problem_boundaries" for the outer boundaries inscribing the problem, "lattice_bounding_box" for the surface bounding the pin lattice, or any other user defined name for regular surfaces.
        - surface_id : unique identifier for the surface
            Rescritions : 0 for the outer boundaries inscribing the problem, -1 for the surface bounding the pin lattice, or any other user defined integer for regular surfaces.
        - center : (x, y) coordinates of the center of the square
        - side_length : the length of the sides of the square
        - corner_radius_of_curvature : if True, the corners of the square will be rounded
        """
        self.name = name
        self.surface_id = surface_id
        self.center = center
        self.side_length = side_length
        self.corner_radius_of_curvature = corner_radius_of_curvature

    def resolve_vertices(self):
        x, y = self.center
        side = self.side_length
        hs = side / 2
        vertices = [
            (x - hs, y - hs),
            (x + hs, y - hs),
            (x + hs, y + hs),
            (x - hs, y + hs)
        ]
        return vertices
    

class Rectangle:
    def __init__(self, name, surface_id, center, width, height, corner_radius_of_curvature=None):
        self.name = name
        self.surface_id = surface_id
        self.center = center
        self.width = width
        self.height = height
        self.corner_radius_of_curvature = corner_radius_of_curvature

    def resolve_vertices(self):
        x, y = self.center
        w_half = self.width / 2
        h_half = self.height / 2
        vertices = [
            (x - w_half, y - h_half),
            (x + w_half, y - h_half),
            (x + w_half, y + h_half),
            (x - w_half, y + h_half)
        ]
        return vertices
    

class Circle:
    def __init__(self, name, surface_id, center, radius):
        self.name = name
        self.surface_id = surface_id
        self.center = center
        self.radius = radius

    def resolve_vertices(self):
        # For a circle, we can return the center and radius as the defining parameters
        return self.center, self.radius
    