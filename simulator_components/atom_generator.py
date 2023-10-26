import pandas as pd
from pandas import DataFrame

from misc_functions import print_debug_key
class AtomGenerator():

    def __init__(self, simulation):
        self.atom_system = simulation.atom_system
        self.prepared_atoms = dict()
        self.settings = simulation.settings

    def prepare(self, atom_name, timestep, position, velocity):
        atom_key = self.atom_system.atom_library[atom_name]['atom_key']
        new_atom = {
        'atom_name': atom_name,
        'atom_key': atom_key,
        'position': position,
        'force': (0, 0, 0),
        'velocity': velocity,
        'status': 'solitary'
        }
        
        if timestep in self.prepared_atoms.keys():
            prepared_atoms_at_timestep = self.prepared_atoms[timestep]
        else:
            prepared_atoms_at_timestep = list()
            
        prepared_atoms_at_timestep.append(new_atom)
        self.prepared_atoms[timestep] = prepared_atoms_at_timestep

    def generate(self, step):
        if step in self.prepared_atoms.keys():
            prepared_atoms = self.prepared_atoms[step]
        else:
            return

        new_atoms = DataFrame(prepared_atoms)
        joined_dataframe = pd.concat([self.atom_system.atoms, new_atoms],
            ignore_index=True)
        self.atom_system.atoms = joined_dataframe
        self.atom_system.update_atom_count()
        print_debug_key(self.settings, 4001)








        
