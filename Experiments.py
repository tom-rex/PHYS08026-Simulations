# -*- coding: utf-8 -*-
"""
File for running Solar System experiments.

This script simulates planetary motion using different numerical methods.
It includes functions for calculating orbital periods, analyzing energy
conservation, and detecting planetary alignments.

Created on Wed Mar 12 13:52:14 2025

@author: Tom Cummings
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from SolarSystem import Celestial_Bodies, Animation

# Define constants
G = 1.18638e-4  # Gravitational constant in appropriate units
M = 332946.0    # Mass of the Sun in solar masses

if __name__ == "__main__":
    # Load simulation parameters from JSON file
    with open('parameters_solar.json') as f:
        parameters_solar = json.load(f)
    
    num_iterations = parameters_solar['num_iterations']  # Number of simulation steps
    timestep = parameters_solar['timestep']  # Time step for integration

    # Initialize planetary data (velocity, mass, radius, color, name)
    planets = [
        [np.array([0, np.sqrt(G * M / body['orbital_radius'])]) if body['orbital_radius'] != 0 else [0, 0],
         body['mass'], body['orbital_radius'], body['colour'], body['name']] for body in parameters_solar['bodies']
    ]

    # Create celestial body system and run simulation
    bodies = Celestial_Bodies(planets, timestep, num_iterations)
    
  #uncomment what integration method you want to use.
    #planets_position = bodies.orbit_direct_euler()
    planets_position = bodies.orbit_euler_cromer()
  #  planets_position = bodies.orbit_beeman()

    # Run animation of planetary motion
    sim = Animation(planets_position, planets)
    sim.run()

def experiment1():
    """
    Experiment 1: Calculate orbital periods of planets.
    """
    bodies.calculate_periods()

def experiment2():
    """
    Experiment 2: Compare energy conservation across different numerical methods.
    """
    # Compute planetary motion using different integration methods
    planets_position = bodies.orbit_euler_cromer()
    KE_cromer, PE_cromer, E_cromer = bodies.calculate_energy()
    
    planets_position = bodies.orbit_direct_euler()
    KE_direct, PE_direct, E_direct = bodies.calculate_energy()
    
    planets_position = bodies.orbit_beehman()
    KE_beehman, PE_beehman, E_beehman = bodies.calculate_energy()
    
    # Generate time array
    T = [timestep * i for i in range(num_iterations)]
    
    # Plot Kinetic Energy
    plt.figure(figsize=(10, 6))
    plt.plot(T, KE_cromer, color="r", label="Euler-Cromer")
    plt.plot(T, KE_direct, color="b", label="Direct Euler")
    plt.plot(T, KE_beehman, color="g", label="Beeman")
    plt.xlabel('Time (years)')
    plt.ylabel('Kinetic Energy (Joules)')
    plt.legend()
    plt.title('Kinetic Energy')
    plt.show()
    
    # Plot Potential Energy
    plt.figure(figsize=(10, 6))
    plt.plot(T, PE_cromer, color="r", label="Euler-Cromer")
    plt.plot(T, PE_direct, color="b", label="Direct Euler")
    plt.plot(T, PE_beehman, color="g", label="Beeman")
    plt.xlabel('Time (years)')
    plt.ylabel('Potential Energy (Joules)')
    plt.legend()
    plt.title('Potential Energy')
    plt.show()
    
    # Plot Total Energy
    plt.figure(figsize=(10, 6))
    plt.plot(T, E_cromer, color="r", label="Euler-Cromer")
    plt.plot(T, E_direct, color="b", label="Direct Euler")
    plt.plot(T, E_beehman, color="g", label="Beeman")
    plt.xlabel('Time (years)')
    plt.ylabel('Total Energy (Joules)')
    plt.legend()
    plt.title('Total Energy')
    plt.show()

def experiment4():
    """
    Experiment 4: Detect planetary alignments within a given angular threshold.
    
    Returns:
        list: A list of times (in years) when planetary alignment occurs.
    """
    planet_positions = bodies.orbit_beehman()[1:7]  # Extract planet positions (excluding the Sun)
    threshold = np.deg2rad(5)  # Alignment threshold in radians
    
    alignment_occurrences = []  # Store times when alignment occurs
    
    for t in range(num_iterations):
        angles = []
        for i in range(5):
            angles.append((np.atan2(planet_positions[i][t][1], planet_positions[i][t][0]) % (np.pi/2)))
        
        mean_angle = sum(angles) / 5
        
        # Check if all angles fall within the threshold of the mean
        if max(angles) < mean_angle + threshold and min(angles) > mean_angle - threshold:
            alignment_occurrences.append(t * timestep)
    
    return alignment_occurrences

        
        