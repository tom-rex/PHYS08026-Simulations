# %%
# %%
import math
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
G = 6.6743 * 10**(-11)

# Mars data
m2 = 6.4185 * 10**23
m2 = 6.4185 * 10**23
r2 = np.array([0, 0])
v2 = np.array([0, 0])

# Phobos data
m1 = 1.06 * 10**16  # Corrected mass (was too low)
r1 = np.array([9.3773 * 10**6, 0])
v1 = np.array([0, np.sqrt(G * m2 / norm(r1))])  # Circular orbit velocity


phobos = [v1, m1, r1]
mars = [v2, m2, r2]


class Celestial_Bodies:
    def __init__(self, planet1, planet2, length_timestep, num_timesteps):
        self.phobos_velocity = np.array(planet1[0], dtype=np.float64)
        self.phobos_mass = planet1[1]
        self.phobos_position = np.array(planet1[2], dtype=np.float64)

        self.mars_velocity = np.array(planet2[0], dtype=np.float64)
        self.mars_mass = planet2[1]
        self.mars_position = np.array(planet2[2], dtype=np.float64)

        self.num_timesteps = num_timesteps
        self.dt = length_timestep

    def orbit(self):
        self.phobos_positions = [self.phobos_position]
        self.phobos_velocities = [self.phobos_velocity]

        self.mars_positions = [self.mars_position]
        self.mars_velocities = [self.mars_velocity]

        for i in range(self.num_timesteps):
            
            kinetic_energy = 0.5*self.mars_mass*norm(self.mars_velocities[-1])**2 +  0.5*self.phobos_mass*norm(self.phobos_velocities[-1])**2
            print(kinetic_energy)
            
            r12 = self.phobos_positions[-1] - self.mars_positions[-1]
            r_mag = norm(r12)

            # Compute gravitational acceleration
            acc_phobos = -G * self.mars_mass / r_mag**3 * r12
            acc_mars = G * self.phobos_mass / r_mag**3 * r12  # Mars also feels force

            # Update velocities using acceleration
            new_phobos_velocity = self.phobos_velocities[-1] + acc_phobos * self.dt
            new_mars_velocity = self.mars_velocities[-1] + acc_mars * self.dt

            # Update positions
            new_phobos_position = self.phobos_positions[-1] + new_phobos_velocity * self.dt
            new_mars_position = self.mars_positions[-1] + new_mars_velocity * self.dt

    
            self.phobos_positions.append(new_phobos_position)
            self.phobos_velocities.append(new_phobos_velocity)

            self.mars_positions.append(new_mars_position)
            self.mars_velocities.append(new_mars_velocity)

        return self.phobos_positions, self.mars_positions


bodies = Celestial_Bodies(phobos, mars, 100, 1000)
phobos_positions, mars_positions = bodies.orbit()

# Convert positions to arrays for plotting
phobos_positions = np.array(phobos_positions) / 1e6  # Convert to Mm
mars_positions = np.array(mars_positions) / 1e6


class Animation:
    def __init__(self):
        
        self.fig, self.ax = plt.subplots()

        self.ax.set_xlim(-12, 12)
        self.ax.set_ylim(-12, 12)
        self.ax.set_xlabel("x (Mm)")
        self.ax.set_ylabel("y (Mm)")
        self.ax.set_aspect("equal")

        # Create planets
        self.phobos_patch = plt.Circle((phobos_positions[0, 0], phobos_positions[0, 1]), 0.5, color="grey", animated=True)
        self.mars_patch = plt.Circle((mars_positions[0, 0], mars_positions[0, 1]), 1, color="r", animated=True)

        self.ax.add_patch(self.phobos_patch)
        self.ax.add_patch(self.mars_patch)

    def animate(self, i):
        self.phobos_patch.center = (phobos_positions[i, 0], phobos_positions[i, 1])
        self.mars_patch.center = (mars_positions[i, 0], mars_positions[i, 1])
        return self.phobos_patch, self.mars_patch

    def run(self):
        self.anim = FuncAnimation(self.fig, self.animate, frames=len(phobos_positions), interval=20, blit=True)
        plt.show()


if __name__ == "__main__":
    
    sim = Animation()
    sim.run()
# %%

# %%