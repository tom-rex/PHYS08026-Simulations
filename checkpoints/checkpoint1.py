# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 18:30:11 2025

@author: Tom Cummings
"""

#Checkpoint 1 task

class Polynomial:
    """A class to represent a 2D vector."""
    def __init__(self, coeffs):
        """constructor setting the coefficients of the polnomial as a list."""
        self.coeffs = coeffs
            
    def order(self):
        """calculate the order and return the result."""
        return len(self.coeffs) - 1
    
    
    def add(self, poly):
        """add a polynomial to self and return the new polynomial"""
        newcoeffs = self.coeffs.copy()
        newpoly = poly.coeffs.copy()
        
        
        while len(newcoeffs) < len(newpoly):
            newcoeffs.append(0)
            
        while len(newpoly) < len(newcoeffs):
            newpoly.append(0)
            
            
        poly_sum = []
        for i in range(0, len(newpoly)):
            poly_sum.append(newcoeffs[i] + newpoly[i])
            

        return Polynomial(poly_sum)

    def derivative(self):
        """Calculate and return the derivative of self"""
        derivative = []
        
        for i in range(1, len(self.coeffs)):
            derivative.append(i*self.coeffs[i])
            
        return Polynomial(derivative)
    
    
    def integral(self, constant):
        """Calculate and return integral of self"""
        integral = [constant]
        
        for i in range(0, len(self.coeffs)):
            integral.append(self.coeffs[i]/(i + 1))
            
        return Polynomial(integral)
    
    
    def __str__(self):
        """string representation in the required format."""
        
        #finding first non-zero term in polynomail
        found = False
        first_nonzero_term = 0
        
        while found == False:
            if self.coeffs[first_nonzero_term] == 0:
                first_nonzero_term += 1
            else:
                found = True
        
        if self.coeffs[0] != 0:
            string = f"{self.coeffs[0]}"
            
            for i in range(1, len(self.coeffs)):
                if self.coeffs[i] != 0:
                    string += f" + {self.coeffs[i]}*x^{i}"
                    
        else:
            string = f"{self.coeffs[first_nonzero_term]}*x^{first_nonzero_term}"
            
            for i in range(first_nonzero_term + 1, len(self.coeffs)):
                if self.coeffs[i] != 0:
                    string += f" + {self.coeffs[i]}*x^{i}"
        return string
                
    
def main():
        """test `Polynomial` class."""
        Pa = Polynomial([2, 0, 4, -1, 0, 6])
        Pb = Polynomial([-1, -3, 0 , 4.5])
        test = Polynomial([0,0,3])
        print(test)
        
        print(f"Pa = {Pa},  Pb = {Pb}")
        print(f"1. The order of Pa(x) is {Pa.order()}")
        print(f"2. Adding Pb(x) to Pa(x) gives {Pa.add(Pb)}")
        print(f"3. The derivative of Pa(x) is {Pa.derivative()} ")
        print(f"4. The integral of the previous expression is {(Pa.derivative()).integral(2)} where c = 2")
        
    
if __name__ == "__main__":
        main()
        

    
    
    
            
        
        
    
