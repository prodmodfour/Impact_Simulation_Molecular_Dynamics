from pandas import Series, DataFrame
import numpy as np
import math

"""
A number of functions used to determine the forces acting on an atom. Each 
function applies to a single atom. Forces between boundary atoms are modelled
as springs whereas forces between bulk atoms and between bulk-boundary atoms are 
Lennard Jones potentials. 
"""

def calculate_spring_constant(epsilon, sigma):
    r = sigma # equillibrium length

    k = ((48 * epsilon) / r **2) * (13 * (sigma / r) - 3.5 * (sigma / r)**6)
    return k

def calculate_x_component(spring_force, atom_variables):
    atoms, index_a, index_b = atom_variables

    position_a = atoms['position'].iloc(index_a)
    position_b = atoms['position'].iloc(index_b)

    a_x, a_y, a_z = position_a
    b_x, b_y, b_z = position_b

    ajacent = a_x - b_x
    opposite = a_y - b_y 

    hypotenuse = math.sqrt((adjacent ^ 2) + (opposite ^ 2))

    cos_theta = adjacent / hypotenuse
    spring_force_x = spring_force * cos_theta
    return spring_force_x

def calculate_y_component(spring_force, atom_variables):
    atoms, index_a, index_b = atom_variables

    position_a = atoms['position'].iloc(index_a)
    position_b = atoms['position'].iloc(index_b)

    a_x, a_y, a_z = position_a
    b_x, b_y, b_z = position_b

    adjacent = a_y - b_y 
    opposite = a_x - b_x

    hypotenuse = math.sqrt((adjacent ^ 2) + (opposite ^ 2))

    cos_theta = adjacent / hypotenuse
    spring_force_y = spring_force * cos_theta
    return spring_force_y

def calculate_z_component(spring_force, atom_variables):
    atoms, index_a, index_b = atom_variables

    position_a = atoms['position'].iloc(index_a)
    position_b = atoms['position'].iloc(index_b)

    a_x, a_y, a_z = position_a
    b_x, b_y, b_z = position_b

    adjacent = a_z - b_z

    x = a_x - b_x
    y = a_y - b_y
    opposite = math.sqrt((x ^ 2) + (y ^ 2))

    hypotenuse = math.sqrt((adjacent ^ 2) + (opposite ^ 2))

    cos_theta = adjacent / hypotenuse
    spring_force_y = spring_force * cos_theta
    return spring_force_z



