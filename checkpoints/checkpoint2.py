#new versin of checkpoint2

#checkpoint 2
import numpy as np
import random


test_lambda = 0.02775
test_N = 20
test_dt = 0.01
#matrix class which has properties such 

class Matrix:
    
    def __init__(self, N):
        self.N = N
        self.atom_matrix = np.ones(shape=(self.N, self.N))
        
    def sum_vals(self):
        
        sum_values = 0
        
        for i in range(self.N):
            for j in range(self.N):     
                
                sum_values += self.atom_matrix[i,j]
        
        return sum_values
    
    def view(self):       
        
        for i in range(0, self.N):
            print(self.atom_matrix[i,:])
        
       # print(str(self.atom_matrix))
        
    def edit_index(self, i, j, val):
        self.atom_matrix[i,j] = val
        
      
    


class Simulate2D(Matrix):
    
    def __init__(self, decay_lambda, N, dt):
        
        self.decay_lambda = decay_lambda
        self.N = N
        self.dt = dt
        super().__init__(N)
        
        
    def decay(self):
                  
        counter = 0
        
        while self.sum_vals() > self.N**2/2 :
            
            for i in range(0,self.N):
                for j in range(0, self.N):
                    
                    RND = random.random()
                    probability = self.decay_lambda * self.dt
                                   
                    if RND <= probability:
                        
                        self.edit_index(i,j, 0)
            counter += 1
                        
        self.view()
        print(f"Intial undecayed nuclei: {self.N*self.N}")
        print(f"Final undecayed nuclei: {self.sum_vals()}")
        print(f"Simulated half life: {counter*self.dt}")
        print(f"Actual half life: {np.log(2)/self.decay_lambda}")

        
        return self
    
    
M = Simulate2D(test_lambda, test_N, test_dt)
M.decay()
    
#M = Matrix(N=20)