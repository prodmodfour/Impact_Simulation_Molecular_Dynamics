import textfile
import os
from misc_functions import print_debug_key



class HistoryWriter:

    def __init__(self, simulation):
        self.simulation = simulation
        self.settings = simulation.settings



    def initialise_trajectory_file(self):
        number_of_atoms = str(self.simulation.atom_system.number_of_atoms)

        filename = "output/trajectories.txt"
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, "w") as file:
            file.write(number_of_atoms)

        with open(filename, 'a') as file:
            file.write('\ntimestep 0') # Title

            self._write_dataframe_to_file(file)

        print_debug_key(self.settings, 3001)
    def _write_dataframe_to_file(self, file):
        atom_dataframe = self.simulation.atom_system.atoms
        for index, row in atom_dataframe.iterrows():
            atom_key = row['atom_key']

            position = row['position']
            x, y, z = position

            line = f"\n{atom_key} {x} {y} {z}"
            file.write(line)

    def append_initial_configuration_settings(self, text):
        textfile('output/initial_configuration.txt', text)

    def write_history(self, atom_dataframe):
        number_of_atoms = str(self.simulation.atom_system.number_of_atoms)
        current_timestep = str(self.simulation.current_timestep)

        filename = "output/trajectories.txt"
        with open(filename, 'a') as file:
            file.write(f"\n{number_of_atoms}")
            file.write(f'\ntimestep {current_timestep}')

            self._write_dataframe_to_file(file)
        print_debug_key(self.settings, 3002)
