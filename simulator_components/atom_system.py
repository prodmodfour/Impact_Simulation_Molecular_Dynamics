import random
import math
from pandas import Series, DataFrame

from misc_functions import check_if_integer, flatten_two_d_array,\
 remove_duplicates, print_debug_key
from unit_cell import UnitCell

class AtomSystem:
    """
    A class used to represent a system of atoms and to perform functions on that
    system.
    """
    def __init__(self, settings):
        self.settings = settings
        
        # position, forces and velocity are tuples. They each contain 3 values
        # and are ordered as x, y, z. E.g. the position tuple contains
        # (x, y, z) position values. Distance is measured in angstroms. Time
        # unit is picoseconds.

        # Status is either 'solitary', 'boundary' or 'bulk'. This represents
        # where in the system the atom is and how the bonds between it and 
        # other atoms is modelled. 
        self.atoms = None
        self.initial_state = None
        self.atom_library = dict()
        self.current_atom_data = dict()
        
        # Physical attributes
        self.width = 0 # Unit Cells
        self.height = 0 # Unit Cells
        self.number_of_layers = 0 # Unit Cells
        self.bulk_range = None

        self.number_of_atoms = 0

        # Energy attributes
        self.temperature = self.settings.required_temperature
        self.kinetic_energy = None
        self.potential_energy = None

    def add_atom_to_library(self, atom_name, atom_data):
        self.atom_library[atom_name] = atom_data

    def set_current_atom_data(self, atom_name):
        self.current_atom_data = self.atom_library[atom_name]

    def set_structure_dimensions(self, dimension_tuple):
        width, height, number_of_layers = dimension_tuple
        check_if_integer(width, height, number_of_layers)
        self.width = width
        self.height = height
        self.number_of_layers = number_of_layers

    def update_atom_count(self):
        self.number_of_atoms = len(self.atoms)

    def reset_energy(self):
        self.kinetic_energy = 0
        self.potential_energy = 0

    def update_energy(self, kinetic_energy, potential_energy):
        self.kinetic_energy = kinetic_energy / (10**10)**2 # To reverse scaling
        self.potential_energy = potential_energy

    def calculate_kinetic_energy(self, index):
        atom_name = self.atoms.at[index, 'atom_name']
        mass = self.atom_library[atom_name]['mass']
        velocity = self.atoms.at[index, 'velocity']
        
        kinetic_energy = 0
        for velocity_component in velocity:
            kinetic_energy += 0.5 * mass * (velocity_component**2)
        if self.settings.debug_keys is True:
            print_debug_key(self.settings, 1009)
        return kinetic_energy

    def calculate_potential_energy(self, index, row):
        return 1.0

    def update_temperature(self):
        boltzmann_constant = 1.380649e-23 # joules per Kelvin
        self.temperature = ((2 * self.kinetic_energy) / 
        (3 * (self.number_of_atoms - 1) * boltzmann_constant))
        print_debug_key(self.settings, 1110)
        print(self.temperature)

    def initialise_atom_structure_pure_metal(self, atom_name):
        """
        Dimension parameters are in units of unit cells. number_of_layers is 
        the number of times that the pattern repeats in the z direction.  Please
        refer to the image 'atom_pattern_pure_metal.png' for more information.
        """
        self.set_current_atom_data(atom_name)
        first_layer = self._create_layer_pure_metal()
        layer_list = self._repeat_layer(first_layer)
        atom_positions = self._convert_to_atom_position_list(layer_list)
        self.populate_dataframe(atom_positions)
        self.initial_state = self.atoms
        self.update_atom_count()

    def _create_layer_pure_metal(self):
        atomic_radius = self.current_atom_data['atomic_radius']
        spacing = self.current_atom_data['spacing']

        # Uses the atom pattern shown in 'atom_pattern_pure_metal.png' 
        layer = list()
        
        # This method creates duplicates of each shared atom. These will be 
        # removed later in the _convert_to_atom_position_list function. 
        top_left_atom = (0, 0, 0)
        for column in range(self.width):
            column = list()

            for cell in range(self.height):
                cell = UnitCell(top_left_atom, atomic_radius, spacing)
                column.append(cell)

                # When building a column of unit cells down, each subsequent
                # cell's top left atom is  the bottom left atom of the previous 
                # cell
                top_left_atom = cell.atoms['bottom_left']
                
                print_debug_key(self.settings, 1001)

            layer.append(column)
            top_cell = layer[0][0]
            # When building a row of columns to the right, each subsequent
            # column has a top cell whose top left atom is the top right atom
            # of the previous cell. 
            top_left_atom = top_cell.atoms['top_right']

        layer = flatten_two_d_array(layer)
        return layer

    def _repeat_layer(self, first_layer):
        layer_list = list()
        layer_list.append(first_layer)

        # All layers are identical. The only difference is that the z
        # coordinate of each atom in each subsequent layer downwards
        # is (2 * atomic_radius) + spacing below the corresponding atom in
        # the preceding layer
        atomic_radius = self.current_atom_data['atomic_radius']
        spacing = self.current_atom_data['spacing']
        delta_z = (2 * atomic_radius) + spacing

        current_layer = first_layer
        for layer in range(1, self.number_of_layers):
            new_layer = list()

            for cell in current_layer:
                top_left_x, top_left_y, top_left_z = cell.atoms['top_left']
                top_left_z -= delta_z
                top_left_atom = (top_left_x, top_left_y, top_left_z)

                new_cell = UnitCell(top_left_atom, atomic_radius, spacing)
                new_layer.append(new_cell)

            layer_list.append(new_layer)
            current_layer = new_layer
            print_debug_key(self.settings, 1002)

        return layer_list

    def _convert_to_atom_position_list(self, layer_list):
        """
        Converts a 2d layer list into a 1d atom position list. Removes duplicate
        atoms.
        """
        cell_list = flatten_two_d_array(layer_list)
        atom_positions = list()
        for cell in cell_list:
            self._add_cell_atom_positions(cell, atom_positions)
            print_debug_key(self.settings, 1003)

        return atom_positions

    def _add_cell_atom_positions(self, cell, atom_positions):
        for label, atom_position in cell.atoms.items():
            if atom_position in atom_positions:
                continue
            elif atom_position not in atom_positions:
                atom_positions.append(atom_position)

    def populate_dataframe(self, atom_positions):
        dataframe_rows = list()
        self.bulk_range = self._determine_bulk_range()

        mass = self.current_atom_data['mass']
        mean = 0
        boltzmann_constant = 1.380649e-23 # joules per Kelvin
        standard_deviation = (boltzmann_constant * self.temperature) / mass
        for atom_position in atom_positions:
            if self.settings.initial_velocity is True:
                velocity = self._randomise_velocity(mean, standard_deviation)
            else:
                velocity = (0, 0, 0)

            new_entry = {
            'atom_name': self.current_atom_data['atom_name'],
            'atom_key': self.current_atom_data['atom_key'],
            'position': atom_position,
            'force': (0, 0, 0),
            'velocity': velocity,
            'status': self._determine_status(atom_position),
            }
            
            dataframe_rows.append(new_entry)
            print_debug_key(self.settings, 1004)
        self.atoms = DataFrame(dataframe_rows)

    def _randomise_velocity(self, mean_velocity, standard_deviation):
        delta_x, delta_y, delta_z = 0, 0, 0
        components = (delta_x, delta_y, delta_z)
        for component in components:
            component =  random.gauss(mean_velocity, standard_deviation)

        velocity = (delta_x, delta_y, delta_z)
        print_debug_key(self.settings, 1006)
        return velocity

    def _determine_bulk_range(self):
        """
        Determines the range of x, y and z values that describe the bulk of the
        model. Any atom that starts in this volume is a bulk atom.
        """
        bulk_range = dict()

        x_boundary_1 = 0 + self.settings.boundary_cutoff
        width = self._calculate_length_angstroms(self.width)
        x_boundary_2 = width - self.settings.boundary_cutoff
        bulk_range['x'] = (x_boundary_1, x_boundary_2)

        y_boundary_1 = 0 - self.settings.boundary_cutoff
        height = self._calculate_length_angstroms(self.height)
        y_boundary_2 = -width + self.settings.boundary_cutoff
        bulk_range['y'] = (y_boundary_1, y_boundary_2)

        z_boundary_1 = 0 - self.settings.boundary_cutoff
        depth = self._calculate_length_angstroms(self.number_of_layers)
        z_boundary_2 = -depth + self.settings.boundary_cutoff
        bulk_range['z'] = (z_boundary_1, z_boundary_2)

        return bulk_range

    def _calculate_length_angstroms(self, number_of_unit_cells):
        atomic_radius = self.current_atom_data['atomic_radius']
        spacing = self.current_atom_data['spacing']

        unit_cell_length = (2 * atomic_radius) + spacing
        length = unit_cell_length * number_of_unit_cells
        return length

    def _determine_status(self, atom_position):
        if self._check_if_in_bulk(atom_position) is True:
            print_debug_key(self.settings, 1011)
            return 'bulk'
        elif self._check_if_in_bulk(atom_position) is False:
            print_debug_key(self.settings, 1022)
            return 'boundary'

    def _check_if_in_bulk(self, atom_position):
        # If it's within the boundaries of each axis, it is in the bulk.
        x, y, z = atom_position
        x_boundary_1, x_boundary_2 = self.bulk_range['x']
        y_boundary_1, y_boundary_2 = self.bulk_range['y']
        z_boundary_1, z_boundary_2 = self.bulk_range['z']

        if (x_boundary_1 <= x <= x_boundary_2) is False:
            return False
        elif (y_boundary_1 >= y >= y_boundary_2) is False:
            return False
        elif (z_boundary_1 >= z >= z_boundary_2) is False:
            return False
        else:
            return True

    def update_atom_forces(self, index):
        self.reset_atom_force(index)
        self._assign_forces(index)
        print_debug_key(self.settings, 1005)

    def reset_atom_force(self, *args):
        for index in args:
            self.atoms.at[index, 'force'] = (0, 0, 0)

    def add_force_to_atom(self, index, force):
        force_x, force_y, force_z = force

        total_force = self.atoms.at[index, 'force']
        total_force_x, total_force_y, total_force_z = total_force

        total_force_x += force_x
        total_force_y += force_y
        total_force_z += force_z
        total_force = (total_force_x, total_force_y, total_force_z)

        self.atoms.at[index, 'force'] = total_force

    def _assign_forces(self, index_a):
        """
        Assign forces on atom a from every other atom b.
        """        
        row_a = self.atoms.iloc[index_a]
        forces = list()
        epsilon = self.current_atom_data['well_depth']
        sigma = self.current_atom_data['van_der_waals']

        for index_b, row_b in self.atoms.iterrows():
            # The atom won't apply forces on itself.
            # We skip calculating forces for atoms that we have already 
            # paired up as we apply the force to both atoms after the first
            # time that they were paired up.
            if index_a >= index_b:
                continue

            status_a = row_a['status']
            status_b = row_b['status']
            if self.settings.boundary_springs is False:
                forces = self._calculate_lennard_jones(row_a, row_b, 
                    epsilon, sigma)
                force_a, force_b = forces

                self.add_force_to_atom(index_a, force_a)
                self.add_force_to_atom(index_b, force_b)
            elif status_a == 'boundary' and status_b == 'boundary':
                forces = self._calculate_spring_force(row_a, row_b)
                force_a, force_b = forces

                self.add_force_to_atom(index_a, force_a)
                self.add_force_to_atom(index_b, force_b)
            else:
                forces = self._calculate_lennard_jones(row_a, row_b, 
                    epsilon, sigma)
                force_a, force_b = forces

                self.add_force_to_atom(index_a, force_a)
                self.add_force_to_atom(index_b, force_b)

    def _calculate_lennard_jones(self, row_a, row_b, epsilon, sigma):
        position_a = row_a['position']
        position_b = row_b['position']
        r = math.dist(position_a, position_b)
        if self.settings.rcutt_off_on and r >= self.settings.rcutt_off:
            force_a, force_b = (0, 0, 0), (0, 0, 0)
            return force_a, force_b

        x_a, y_a, z_a = position_a
        x_b, y_b, z_b = position_b

        r_x = x_a - x_b
        r_y = y_a - y_b
        r_z = z_a - z_b
        components = (r_x, r_y, r_z)

        force_components_a = list()
        for r_comp in components:
            force = ((24 * epsilon) / r**2) * ((sigma / r)**12 - (sigma / r)**6) * r_comp
            force_components_a.append(force)
        force_a = tuple(force_components_a)

        force_components_b = list()
        for component in force_components_a:
            component = -component
            force_components_b.append(component)
        force_b = tuple(force_components_b)
        
        print_debug_key(self.settings, 1007)
        return force_a, force_b

    def _calculate_spring_force(self, row_a, row_b):
        position_a = row_a['position']
        position_b = row_b['position']
        r = math.dist(position_a, position_b)

        if self.settings.rcutt_off_on and r >= self.settings.rcutt_off:
            force_a, force_b = (0, 0, 0), (0, 0, 0)
            return force_a, force_b

        atom_name = row_a['atom_name']
        r_equillibrium = self.atom_library[atom_name]['van_der_waals']
        spring_constant = self.atom_library[atom_name]['spring_constant']

        x_a, y_a, z_a = position_a
        x_b, y_b, z_b = position_b

        r_x = x_a - x_b
        r_y = y_a - y_b
        r_z = z_a - z_b
        components = (r_x, r_y, r_z)

        force_components_a = list()
        for r_comp in components:
            force = (-spring_constant * (r - r_equillibrium)) * r_comp
            force_components_a.append(force)
        force_a = tuple(force_components_a)

        force_components_b = list()
        for component in force_components_a:
            component = -component
            force_components_b.append(component)
        force_b = tuple(force_components_b)
        
        print_debug_key(self.settings, 1008)
        return force_a, force_b
