#################################################################
# Name:     PMGrid.py                                           #
# Authors:  Yuan Qi Ni                                          #
#           Michael Battaglia                                   #
# Date:     November 23, 2018                                   #
# Function: Create grid object to be used for mesh in a         #
#           particle mesh N body simulation.	                #
#################################################################

#essential imports
import numpy as np
from gradient import Grad
import matplotlib.pyplot as plt

#class: grid which can hold massive objects
class Grid:
    """array of points as weighted grid"""
    def __init__(self, array, rx, ry, D):
        self.array = array
        self.D = D #D is the length of grid spacings
        self.r = np.array([rx, ry])
    def resetGrid(self):
        #resets grid to zeros
        self.array = np.zeros(self.array.shape)
    def insertBody(self, body):
        #inserts a body's mass into a gridpoint
        pos = np.rint((body.r-self.r)/self.D).astype(int)
        self.array[pos[0],pos[1]] += body.m

    def evalForce(self):
        #density on grid points
        density = self.array/self.D**2
        #size of grid
        M = self.array.shape[0]
        
        #Fourier exponential 
        W = np.exp(2*np.pi*1j/M)

        #Fourier transform of density
        density_fft = (1.0/M)*np.fft.rfft2(density)
        
        #Fourier transform of kernel
        denom = np.zeros(density_fft.shape).astype(complex)
        for i in range(density_fft.shape[0]):
            for j in range(density_fft.shape[1]):
                #kernel value at each frequency
                denom[i][j] = -(W**i+W**-i+W**j+W**-j-4)/self.D**2

        #Fourier transform of potential
        potential_fft = density_fft/denom
        potential_fft[0][0] = 0.0 + 0.0j #have to set this to not get NaNs from dividing by zero
        #potential over grid
        potential_grid = np.fft.irfft2(potential_fft)

        '''
        #Testing scheme to plot potential, ensuring Poisson equation was
        #solved properly. Not part of main code. Run one at a time.
        vmax = np.amax(potential_grid)
        vmin = np.amin(potential_grid)*0.5
        plt.imshow(potential_grid, cmap='Greys', vmax=vmax, vmin=vmin, interpolation='nearest')
        plt.show()
        '''

        #force gradient of potential over grid
        force_x, force_y = Grad(potential_grid, self.D)
        force_grid = np.transpose(np.array([force_x,force_y]), (1,2,0))
        #update forces on grid points
        self.forces = force_grid
    def forceOn(self, body):
        #round positions on grid to nearest int
        pos = np.rint((body.r-self.r)/self.D).astype(int)
        return self.forces[pos[0],pos[1]]

    def plot(self, u):
        #plot grid, and number of masses in each grid
        num = (self.array/u).astype(int)
        for i in range(num.shape[0]):
            for j in range(num.shape[1]):
                n = num[i][j]
                x = i*D + self.r[0]
                y = j*D + self.r[1]
                plt.text(x, y, n)
