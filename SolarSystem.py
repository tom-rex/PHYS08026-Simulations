# -*- coding: utf-8 -*-
"""
Solar System Simulation using Euler and Beeman's Methods

Created on Wed Mar 12 13:52:14 2025

@author: Tom Cummings
"""

import json
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
G = 1.18638e-4  # Gravitational constant in appropriate units
M = 332946.0    # Mass of the Sun in solar masses


class Celestial_Bodies:
    """
    Class representing a system of celestial bodies, simulating gravitational interactions
    between planets and the Sun using numerical integration methods (Euler and Beeman's).

    Attributes:
        num_planets (int): Number of planets in the system.
        num_timesteps (int): Total number of simulation timesteps.
        dt (float): Length of each timestep.
        velocities (np.array): Initial velocities for each planet.
        positions (np.array): Initial positions for each planet.
        masses (np.array): Masses of each planet.
        names (list): Names of each planet.
    """

    def __init__(self, planets, length_timestep, num_timesteps):
        """
        Initialize the simulation with the given planets, timestep, and number of timesteps.

        Parameters:
            planets (list): List of planets with initial conditions (velocity, mass, radius, color, name).
            length_timestep (float): Length of each simulation timestep.
            num_timesteps (int): Number of timesteps in the simulation.
        """
        self.num_planets = len(planets)
        self.num_timesteps = num_timesteps
        self.dt = length_timestep

        # Initialize arrays to store velocities, positions, and masses
        self.velocities = np.zeros((self.num_planets, 2))
        self.positions = np.zeros((self.num_planets, 2))
        self.masses = np.zeros(self.num_planets)
        self.names = []

        # Populate initial conditions for each planet
        for i in range(self.num_planets):
            self.velocities[i] = np.array(planets[i][0])  # Initial velocities (vx, vy)
            self.masses[i] = float(planets[i][1])  # Mass of the planet
            self.positions[i] = np.array([planets[i][2], 0])  # Initial position along x-axis
            self.names.append(planets[i][4])  # Name of the planet

    def calculate_acceleration(self, i, time):
        """
        Calculate gravitational acceleration acting on planet 'i' due to all other planets.

        Parameters:
            i (int): Index of the planet to calculate acceleration for.
            time (int): Current timestep index.

        Returns:
            np.array: Acceleration vector for planet 'i' at the specified time.
        """
        if time == -1:
            # Special case for initial acceleration calculation
            return self.calculate_acceleration(i, time + 1)
        else:
            acceleration = np.array([0.0, 0.0])

            # Loop through all planets to sum gravitational forces
            for j in range(self.num_planets):
                if j != i:
                    # Calculate vector from planet j to planet i
                    rij = self.planet_positions[i][time] - self.planet_positions[j][time]
                    r_mag = norm(rij)  # Magnitude of separation vector
                    acceleration += -G * self.masses[j] / r_mag**3 * rij  # Gravitational acceleration
            
            return acceleration

    def calculate_energy(self):
        """
        Calculate and print the kinetic energy of the system at every 100 timesteps.

        Returns:
            list: Kinetic energy values at each timestep.
        """
        kinetic_energy = []

        # Loop through each timestep to compute kinetic energy
        for t in range(self.num_timesteps):
            KE = 0
            for i in range(self.num_planets):
                KE += 0.5 * self.masses[i] * norm(self.planet_velocities[i][t])**2
            kinetic_energy.append(KE)
        
        potential_energy = []
        for t in range(self.num_timesteps):
            PE = 0
            for i in range(self.num_planets):
                for j in range(self.num_planets):
                    if i != j:
                        PE +=  -0.5 * G * self.masses[i]*self.masses[j] / norm(self.planet_positions[i][t] - self.planet_positions[j][t])
            potential_energy.append(PE)
        
        
        energy = [kinetic_energy[i] + potential_energy[i] for i in range(self.num_timesteps)]
        # Print energy every 100 timesteps
        for t in np.arange(0, self.num_timesteps, 100):
            print(f"KE_{t} = {kinetic_energy[t]}, PE_{t} = {potential_energy[t]}, Energy_{t} = {energy[t]}")

        return kinetic_energy, potential_energy, energy
    
    def orbit_euler_cromer(self):
        """
        Simulate orbits using Euler-Cromer method:
        v(t+dt) = v(t) + a(t)*dt
        r(t+dt) = r(t) + v(t+dt)*dt
        """
        # Initialize lists to store positions and velocities over time
        self.planet_positions = [[self.positions[i].copy()] for i in range(self.num_planets)]
        self.planet_velocities = [[self.velocities[i].copy()] for i in range(self.num_planets)]
    
        for t in range(self.num_timesteps):
            for i in range(self.num_planets):
                # Euler-Cromer update: first velocity then position
                new_velocity_i = self.planet_velocities[i][-1] + self.calculate_acceleration(i, t) * self.dt
                new_position_i = self.planet_positions[i][-1] + new_velocity_i * self.dt

                self.planet_positions[i].append(new_position_i)
                self.planet_velocities[i].append(new_velocity_i)
                
        return np.array(self.planet_positions)
    

    def orbit_beeman(self):
         """
         Simulate orbits using Beeman's method:
         r(t+dt) = r(t) + v(t)*dt + (1/6)*(4*a(t) - a(t-1))*dt^2
         v(t+dt) = v(t) + (1/6)*(2*a(t+dt) + 5*a(t) - a(t-1))*dt
         """
         # Initialize lists to store positions and velocities over time
         self.planet_positions = [[self.positions[i].copy()] for i in range(self.num_planets)]
         self.planet_velocities = [[self.velocities[i].copy()] for i in range(self.num_planets)]

         # Calculate initial accelerations
         accelerations = [self.calculate_acceleration(i, 0) for i in range(self.num_planets)]
         
         # Calculate previous accelerations (needed for Beeman's method)
         # For the first step, we can use a simple Euler step backwards
         prev_accelerations = []
         for i in range(self.num_planets):
             prev_pos = self.positions[i] - self.velocities[i] * self.dt
             # Store current positions temporarily
             current_pos = self.planet_positions[i][0].copy()
             self.planet_positions[i][0] = prev_pos
             prev_acc = self.calculate_acceleration(i, -1)
             # Restore current positions
             self.planet_positions[i][0] = current_pos
             prev_accelerations.append(prev_acc)

         # Main simulation loop
         for t in range(self.num_timesteps):
             new_accelerations = []

             # Predictor step: predict new positions
             for i in range(self.num_planets):
                 new_position_i = (
                     self.planet_positions[i][-1] +
                     self.planet_velocities[i][-1] * self.dt +
                     (1/6) * (4 * accelerations[i] - prev_accelerations[i]) * self.dt**2
                 )
                 self.planet_positions[i].append(new_position_i)

             # Calculate new accelerations using updated positions
             for i in range(self.num_planets):
                 new_accelerations.append(self.calculate_acceleration(i, t + 1))

             # Corrector step: update velocities using new accelerations
             for i in range(self.num_planets):
                 new_velocity_i = (
                     self.planet_velocities[i][-1] +
                     (1/6) * (2 * new_accelerations[i] + 5 * accelerations[i] - prev_accelerations[i]) * self.dt
                 )
                 self.planet_velocities[i].append(new_velocity_i)

             # Update accelerations for the next step
             prev_accelerations = accelerations.copy()
             accelerations = new_accelerations.copy()

         return np.array(self.planet_positions)
    
    def orbit_direct_euler(self):
        """
        Simulate orbits using Direct Euler method:
        r(t+dt) = r(t) + v(t)*dt
        v(t+dt) = v(t) + a(t)*dt
        """
        # Initialize lists to store positions and velocities over time
        self.planet_positions = [[self.positions[i].copy()] for i in range(self.num_planets)]
        self.planet_velocities = [[self.velocities[i].copy()] for i in range(self.num_planets)]
        
        for t in range(self.num_timesteps):
            for i in range(self.num_planets):
                # Direct Euler update: first position then velocity
                new_position_i = self.planet_positions[i][-1] + self.planet_velocities[i][-1] * self.dt
                new_velocity_i = self.planet_velocities[i][-1] + self.calculate_acceleration(i, t) * self.dt
                
                self.planet_positions[i].append(new_position_i)
                self.planet_velocities[i].append(new_velocity_i)
        
        return np.array(self.planet_positions)
    
    
    def calculate_periods(self):  
        """ Calculate orbital periods for each planet."""
        orbital_periods = []
        for i in range(self.num_planets):
            try:
                 t = 0
                 # since we can't calucalte the orbit of the sun around the sun
                 if (self.planet_positions[i][0][0]) == float(0) and (self.planet_positions[i][0][1] == float(0)):
                     t = 0
                     sun_positions = self.planet_positions[i]
                 else:
                     while (self.planet_positions[i][t][1] < sun_positions[t][1] and self.planet_positions[i][t + 1][1] > sun_positions[t+1][1]) == False:
                         t+=1
                 orbital_periods.append(t*self.dt)
                 print(self.names[i], "has period ", orbital_periods[i])
            except IndexError:
                print("Integration time is not long enough to calculate any further orbits")
                break    

# Animation class to visualize the orbits
class Animation(Celestial_Bodies):
    """
    Class for visualizing planet orbits in 2D space using matplotlib animations.
    """

    def __init__(self, planets_position, planets):
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim(-10, 10)
        self.ax.set_ylim(-10, 10)
        self.ax.set_xlabel("x (Mm)")
        self.ax.set_ylabel("y (Mm)")
        self.ax.set_aspect("equal")
        
        self.planet_positions = planets_position
        self.planets = planets
        
        # Create planet patches for animation
        self.planet_patches = [
            plt.Circle((planets_position[i][0, 0], planets_position[i][0, 1]), 0.1, color=planets[i][3], animated=True)
            for i in range(len(self.planets)) # note that num_planets = len(planets_position)
        ]
        for patch in self.planet_patches:
            self.ax.add_patch(patch)
            
    def animate(self, i):
        """
        Update planet positions for each animation frame.
        """
        for n in range(len(self.planets)):
            self.planet_patches[n].center = (self.planet_positions[n][i, 0],self.planet_positions[n][i, 1])
        return self.planet_patches

    def run(self):
        """
        Run the matplotlib animation.
        """
        self.anim = FuncAnimation(self.fig, self.animate, frames=len(self.planet_positions[0]), interval=20, blit=True)
        plt.show()



