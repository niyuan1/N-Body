#################################################################
# Name:     Grid.py                                             #
# Authors:  Michael Battaglia                                   #
# Course:   Phy407                                              #
# Date:     December 13, 2016                                   #
# Function: Create grid object to be used for mesh in a particle#
#			mesh-style N body simulation.	                    #
#################################################################

#essential imports
import numpy as np
from gradient import Grad #old code from lab 3
import matplotlib.pyplot as plt

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
		#from Body.py
		pos = np.rint((body.r-self.r)/self.D).astype(int)
		self.array[pos[0],pos[1]] += body.m

	def evalForce(self):
		#Use density at points to compute gravitational potentials
		#by solving Poisson's equation with FFTs. 
		#From these potentials, use fft to go to Fourier space so
		#it's easier to apply Green's function, then transform
		#back to our origin space.
		density = self.array/self.D**2
		M = self.array.shape[0]

		W = np.exp(2*np.pi*1j/M) #to write denom more compactly

		density_fft = (1.0/M)*np.fft.rfft2(density) #go to Fourier space
		denom = np.zeros(density_fft.shape).astype(complex) #was converting to reals before, had to specify complex
		for i in range(density_fft.shape[0]):
			for j in range(density_fft.shape[1]):
				#equation from Fourier transform reference (Gonsalves, 2011) in paper
				denom[i][j] = -(W**i+W**-i+W**j+W**-j-4)/self.D**2

		potential_fft = density_fft/denom
		potential_fft[0][0] = 0.0 + 0.0j #have to set this to not get NaNs from dividing by zero
		potential_grid = np.fft.irfft2(potential_fft) #go back to original space


		'''
		#Testing scheme to plot potential, ensuring Poisson equation was
		#solved properly. Not part of main code. Run one at a time.
		vmax = np.amax(potential_grid)
		vmin = np.amin(potential_grid)*0.5
		plt.imshow(potential_grid, cmap='Greys', vmax=vmax, vmin=vmin, interpolation='nearest')
		plt.show()
		'''


		#find the force at the grid points from the potentials by taking
		#the negative of the gradient. Gradient code from lab 3.
		force_x, force_y = Grad(potential_grid, self.D)
		force_grid = np.transpose(np.array([force_x,force_y]), (1,2,0))

		self.forces = force_grid

	def forceOn(self, body):
		pos = np.rint((body.r-self.r)/self.D).astype(int)
		return self.forces[pos[0],pos[1]]

	def plot(self, u):
		#Used to create some plots included but not the main animation.
		num = (self.array/u).astype(int)
		for i in range(num.shape[0]):
			for j in range(num.shape[1]):
				n = num[i][j]
				x = i*D + self.r[0]
				y = j*D + self.r[1]
				plt.text(x, y, n)
