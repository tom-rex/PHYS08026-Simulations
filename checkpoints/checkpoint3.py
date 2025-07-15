#checkpoint 3
import numpy as np
import matplotlib.pyplot as plt

class Mandelbrot:
        
    def convergence(self, real, imaginary):
       
        C = complex(real, imaginary)
        z = 0      
        converges = True
        N = 0
        
        while converges == True and N < 256:
            
            if abs(z) > 2:
                converges = False
                break
            
            else:
                N += 1
                z = z**2 + C
            
        return(256 - N)
    
       
    def plot(self, x_lower, x_upper, y_lower, y_upper):
        
        vectorized_convergence = np.vectorize(self.convergence)

        X = np.linspace(-2.025, 0.6, 512)
        Y = np.linspace(-1.125, 1.125, 512)

        XX, YY = np.meshgrid(X,Y)

        # calculate z values using vectorized function
        Z = vectorized_convergence(XX, YY)



        plt.imshow(Z, extent=(X.min(), X.max(), Y.min(), Y.max()), 
                   vmin=0, vmax=256, cmap='bone', aspect='auto')

       # frame1.axes.get_xaxis().set_visible(False)
        #frame1.axes.get_yaxis().set_visible(False)

        plt.show()

if __name__ == "__main__":  
    x_lower = -2 ; x_upper = 0.6
    y_lower = -1.125 ; y_upper = 1.125
    Mandelbrot().plot(x_lower, x_upper, y_lower, y_upper)


             
            
            
        


