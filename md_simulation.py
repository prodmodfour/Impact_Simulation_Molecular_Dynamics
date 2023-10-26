import math
import random

from simulator_components.atom_system import AtomSystem
from simulator_components.atom_generator import AtomGenerator
from simulator_components.terminal import Terminal
from simulator_components.history_writer import HistoryWriter

from physics.velocity_verlet import calculate_next_half_step_velocity, \
calculate_next_position
from misc_functions import print_debug_key

class MDSimulation():
    """
    Used to manage a molecular dynamics simulation. A simulation is composed of 
    an atomic system, an atom generator, a terminal and a logger.
    """
    def __init__(self, settings):
        self.settings = settings
        self.atom_system = AtomSystem(self.settings)
        self.atom_generator = AtomGenerator(self)
        self.terminal = Terminal(self.settings)
        self.writer = HistoryWriter(self)

        self.current_timestep = 0
        self.simulation_time = 0
        self.history_interval = 10

        self.boltzmann_constant = self.settings.boltzmann_constant
        self.required_temperature = self.settings.required_temperature
        
        # This is for atoms at the boundary
        self.standard_deviation_velocity_copper = None

    def initialise_simulation(self, atom_name, dimension_tuple):
        self.atom_system.set_dimensions(dimension_tuple)
        self.atom_system.initialise_atoms_pure_metal(atom_name)

    def initialise_standard_deviation_velocity_copper(self, mass):
        self.standard_deviation_velocity_copper = (
            (self.boltzmann_constant * self.required_temperature) / mass)

    def run(self, number_timesteps):
        for step in range(number_timesteps):
            self.current_timestep += 1
            self._simulate_step(step)
            print(step)
            print_debug_key(self.settings, 2000)

    def _simulate_step(self, step):
        self.simulation_time += self.settings.timestep_size
        self.atom_generator.generate(step)

        kinetic_energy = 0
        potential_energy = 0
        for index, row in self.atom_system.atoms.iterrows():
            self._update_atom(index)

            if self.atom_system.atoms.at[index, 'status'] == 'boundary':
                self._thermostat(index)

            kinetic_energy += self.atom_system.calculate_kinetic_energy(index)
            potential_energy += self.atom_system.calculate_potential_energy(
                index, row)
        self.atom_system.update_energy(kinetic_energy, potential_energy)
        self.atom_system.update_temperature()

        if self.current_timestep % self.history_interval == 0:
            self.writer.write_history(self.atom_system.atoms)


    def _update_atom(self, index):
        self.atom_system.update_atom_forces(index)

        force = self.atom_system.atoms.at[index, 'force']
        velocity = self.atom_system.atoms.at[index, 'velocity']
        position = self.atom_system.atoms.at[index, 'position']
        atom_name = self.atom_system.atoms.at[index, 'atom_name']
        self.atom_system.set_current_atom_data(atom_name)
        mass = self.atom_system.current_atom_data['mass']
        delta_time = self.settings.timestep_size

        # Velocity Verlet Algorithm
        velocity = calculate_next_half_step_velocity(force, velocity,
            mass, delta_time)
        position = calculate_next_position(position, velocity, delta_time)
        self.atom_system.atoms.at[index, 'position'] = position
        
        self.atom_system.update_atom_forces(index)
        force = self.atom_system.atoms.at[index, 'force']
        velocity = calculate_next_half_step_velocity(force, velocity, mass,
            delta_time)
        self.atom_system.atoms.at[index, 'velocity'] = velocity

    def _thermostat(self, index):
        mean_velocity = 0
        velocity_x, velocity_y, velocity_z = 0, 0, 0
        velocity_components = [velocity_x, velocity_y, velocity_z]
        for velocity_component in velocity_components:
            velocity_component = random.gauss(mean_velocity, 
                self.standard_deviation_velocity_copper)
        
        velocity_components = tuple(velocity_components)
        self.atom_system.atoms.at[index, 'velocity'] = velocity_components
        print_debug_key(self.settings, 2100)
