import math

from md_simulation import MDSimulation
from simulator_components.settings import Settings
from physics.forces import calculate_spring_constant

# Simulation Settings
settings = Settings()
simulation = MDSimulation(settings)

simulation.settings.initial_velocity = False
simulation.settings.rcutt_off_on = True
simulation.settings.boundary_springs = True
simulation.settings.thermostatting = False 

simulation.settings.debug_keys = True

# Atom Data
# We will be firing an argon atom at a block of copper.
mass = 6.63352088e-26 * (10**10) # Convert to same scale as Angstrom
argon_data = {
    'atom_name': 'argon',
    'atom_key': 'Ar',
    'atomic_radius': 0.98, # Angstroms
    'well_depth': 0.010333179484427387, #eV
    'van_der_waals': 1.88, # Angstroms
    'mass': mass # 10^10 kg
}
simulation.atom_system.add_atom_to_library('argon', argon_data)

mass = 1.055e-25 * (10**10) # Convert to same scale as Angstrom
well_depth = 0.4802 # eV
van_der_waals = 2.285 # Angstroms
spring_constant = calculate_spring_constant(well_depth, van_der_waals)

copper_data = {
    'atom_name': 'copper',
    'atom_key': 'Cu',
    'atomic_radius': 1.28, # Angstroms
    'spacing': 3.61, # Angstroms
    'spring_constant': spring_constant,
    'well_depth': well_depth, # ev
    'van_der_waals': van_der_waals, # Angstroms
    'mass': mass # 10^10 kg
}
simulation.atom_system.add_atom_to_library('copper', copper_data)
simulation.initialise_standard_deviation_velocity_copper(mass)

# Initialise Copper Block
width = 20
height = 20
number_of_layers = 25
dimensions = (width, height, number_of_layers)
simulation.atom_system.set_structure_dimensions(dimensions)
simulation.atom_system.initialise_atom_structure_pure_metal('copper')
print(simulation.atom_system.number_of_atoms)

# Setting rcutoff
rcutoff = 2.5 * copper_data['van_der_waals']
simulation.settings.rcutoff = rcutoff

# Setting up trajectory file
simulation.writer.initialise_trajectory_file()

# Set up particles to be generated and fired at specific time steps
timestep = 4
position = (width / 2, height / 2, 10)

argon_mass = simulation.atom_system.atom_library['argon']['mass']
kinetic_energy = 1.60218e-16 # 1k eV in joules
speed = math.sqrt((2 * kinetic_energy) / argon_mass)
speed = speed * 10**10 # Convert to Angstroms / second
velocity = (0, 0, -speed)

simulation.atom_generator.prepare('argon', timestep, position, velocity)

# Run the simulation
number_timesteps = 10
simulation.run(number_timesteps)