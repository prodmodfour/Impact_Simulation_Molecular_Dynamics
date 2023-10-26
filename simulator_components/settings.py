
class Settings:
    def __init__(self):
        """Initialise default settings"""
        self.required_temperature = 300 # Kelvin
        self.timestep_size = 1e-5 # 1 Femtosecond in 10^-10 scale
        self.boltzmann_constant = 1.380649e-23 # joules per Kelvin
        

        # Distance in Angstroms that we consider the boundary of the block.
        # It should be noted that we are not using Period Boundary Conditions.
        # Boundary refers to the outermost section of atoms. 
        self.boundary_cutoff = 2 # Angstroms
        self.rcutt_off = 10 # Angstroms

        # Molecular dynamics settings
        self.initial_velocity = False
        self.rcutt_off_on = False 
        self.boundary_springs = False
        self.thermostatting = False 

        self.debug_keys = False