

class UnitCell:
    """
    Contains methods for the creation and storage of one unit cell of an atom 
    pattern. A unit cell an indivisible portion of a repeating pattern of atoms.
    """

    def __init__(self, top_left_atom, atomic_radius, spacing):
        """
        Atom variables, such at top_left_atom, are tuples with x, y z 
        coordinates. Values are floats.
        """
        # Please refer to the image 'atom_pattern_pure_metal.png' for more info. 
        self.atomic_radius = atomic_radius
        self.spacing = spacing

        self.atoms = {
        'top_left': top_left_atom,
        'bottom_left': None,
        'top_right': None,
        'bottom_right': None,
        'center': None
        }

        self.atoms['bottom_left'] = self._calculate_bottom_left_position()
        self.atoms['top_right'] = self._calculate_top_right_position()
        self.atoms['bottom_right'] = self._calculate_bottom_right_position()
        self.atoms['center'] = self._calculate_center_position()

    def _calculate_bottom_left_position(self):
        top_left_x, top_left_y, top_left_z = self.atoms['top_left']

        bottom_left_x = top_left_x
        bottom_left_y = top_left_y - (self.atomic_radius * 2) - self.spacing
        bottom_left_z = top_left_z

        bottom_left_atom = (bottom_left_x, bottom_left_y, bottom_left_z)
        return bottom_left_atom

    def _calculate_top_right_position(self):
        top_left_x, top_left_y, top_left_z = self.atoms['top_left']

        top_right_x = top_left_x + (self.atomic_radius * 2) + self.spacing
        top_right_y = top_left_y
        top_right_z = top_left_z

        top_right_atom = (top_right_x, top_right_y, top_right_z)
        return top_right_atom

    def _calculate_bottom_right_position(self):
        top_right_x, top_right_y, top_right_z = self.atoms['top_right']

        bottom_right_x = top_right_x
        bottom_right_y = top_right_y - (self.atomic_radius * 2) - self.spacing
        bottom_right_z = top_right_z

        bottom_right_atom = (bottom_right_x, bottom_right_y, bottom_right_z)
        return bottom_right_atom

    def _calculate_center_position(self):
        top_left_x, top_left_y, top_left_z = self.atoms['top_left']
        bottom_left_x, bottom_left_y, bottom_left_z = self.atoms['bottom_left']
        top_right_x, top_right_y, top_right_z = self.atoms['top_right']

        midpoint_x = (top_left_x + top_right_x) / 2
        midpoint_y = (top_left_y + bottom_left_y) / 2
        # Z coordinates for atoms in the same unit cell are always the same
        midpoint_z = top_left_z

        center_atom = (midpoint_x, midpoint_y, midpoint_z)
        return center_atom
        








