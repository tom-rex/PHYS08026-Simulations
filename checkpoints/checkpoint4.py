# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:41:09 2025

@author: Tom Cummings
"""

import numpy as np
import random
import matplotlib.pyplot as plt
        
class Traffic:
    
    def __init__(self, road_size, car_density, num_iterations):
        self.road = np.zeros((num_iterations, road_size), dtype=int)
        for i in range(road_size):
            if random.random() < car_density:
                self.road[0, i] = 1  # Initialize first row with cars
            
        self.iterations = num_iterations
        self.road_length = road_size
        self.num_cars = sum(self.road[0,:])
        self.density = self.num_cars/self.road_length
    
    def simulate(self):
        
        average_speeds = np.zeros(self.iterations)
        
        print("Initial state:\n", self.road[0] ,"speed = 0")
        
        for i in range(self.iterations - 1):  # Iterate over time steps
            new_road = self.road[i].copy()  # Copy current state
            row_speed = 0
            
            for j in range(self.road_length):  # Iterate over road positions
            
                
                next_pos = (j + 1) % self.road_length  # Circular road movement
                prev_pos = (j - 1) % self.road_length

                if self.road[i, j] == 1:
                    if self.road[i, next_pos] == 0:  # Move car forward if space is empty
                        new_road[next_pos] = 1
                        new_road[j] = 0
                        row_speed += 1 #when a car moves forward row speed gains one
                elif self.road[i,j] == 0:
                    if self.road[i, prev_pos] == 1:
                        new_road[next_pos] = 1
                       # row_speed+=1
            #print(row_speed)
            self.road[i + 1] = new_road  # Update the next row with the new state
                   
            if self.num_cars > 0:
                #average_speeds[i] = row_speed/self.num_cars
                average_speeds[i] = row_speed/self.density#putting the speeds into an array
            else:
                average_speeds[i] = 0
                
            print(f"Iterationnn {i+1}:\n", self.road[i + 1], "speed = ", row_speed)  # Print each iteration
             
        overall_average_speed = sum(average_speeds)/self.iterations
        steady_speed = average_speeds[-2]
        print(average_speeds)
        return steady_speed
    
    def steadyDensity(self):
        steadys = []
        densities = np.linspace(0, 1, 1000)
        for i in range(100):
            traffic = Traffic(self.road_length, densities[i], self.iterations)
            steadys.append(traffic.simulate())
            
            
        plt.plot(densities, steadys, marker = "o")
        plt.xlabel('density')
        plt.show()
        
            
            
    
traffic = Traffic(1000,0.3,1000)
result = traffic.simulate()
plot = traffic.steadyDensity()
    

