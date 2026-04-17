### class defining a cell in the geometry :
# assign a name, unique id and list of surfaces with their logical operators


class Cell:
    def __init__(self, name, cell_id, surfaces_to_operators):
        """
        Define a cell with :
        - name : user friendly name for the cell
        - cell_id : unique identifier for the cell
        - surfaces_to_operators : dictionary mapping surface IDs to their logical operators
          (e.g., {1: "+", 2: "-", 3: "#"}), where '+' indicates the outside of the surface, '-' indicates the inside, and '#' indicates a complement.
          
        """
        self.name = name
        self.cell_id = cell_id
        self.surfaces_to_operators = surfaces_to_operators

    def _validate_operators(self):
        valid_operators = {"+", "-", "#"}
        for operator in self.surfaces_to_operators.values():
            if operator not in valid_operators:
                raise ValueError(f"Invalid operator '{operator}' in cell '{self.name}'. "
                                 f"Valid operators are: {valid_operators}")