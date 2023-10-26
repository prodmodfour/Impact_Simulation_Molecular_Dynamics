def calculate_next_half_step_velocity(force, velocity, mass, delta_time):
    delta_x, delta_y, delta_z = velocity
    force_x, force_y, force_z = force
    components = [(delta_x, force_x), (delta_y, force_y), (delta_z, force_z)]

    next_velocity_components = list()
    for component in components:
        velocity, force = component

        velocity_component = velocity + ((force * delta_time) / (2 * mass))
        next_velocity_components.append(velocity_component)

    next_half_step_velocity = tuple(next_velocity_components)
    return next_half_step_velocity

def calculate_next_position(position, velocity, delta_time):
    x, y, z = position
    delta_x, delta_y, delta_z = velocity
    components = [(x, delta_x), (y, delta_y), (z, delta_z)]

    next_position_components = list()
    for component in components:
        position, velocity = component
        next_position = position + (velocity * delta_time)
        next_position_components.append(next_position)

    next_position = tuple(next_position_components)
    return next_position

