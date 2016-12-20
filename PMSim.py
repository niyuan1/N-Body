#################################################################
# Name:     ParticleMesh.py                                     #
# Authors:  Michael Battaglia                                   #
# Course:   Phy407                                              #
# Date:     December 15, 2016                                   #
# Function: Model the dynamics of an N-body system with			#
#			gravitational potentials using a particle mesh 		#
#			method.	Doesn't include testing aglorithms.			#   
#################################################################

#essential imports
import numpy as np
from PMGrid import Grid 
from MCgalaxy import generateGalaxy
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#Constants: Milky Way parameters; decided on with Chris Ni
r0 = 3 #kiloparsecs, scale length of galaxy
m0 = 50.0 #10^9 solar masses, mass of galaxy
#generate 1000 masses in 15kpc box
N = 1000 #number of bodies
L = 15.0 #Number of grid divisions, decided on after trying several divisions. 

#Create bodies data from MCgalaxy.py
bodies = generateGalaxy(r0, m0, N, L)

#grid resolution and initializing grid
D = L/np.sqrt(N) #kpc grid spacing, very fine based on number of bodies in our grid
init = np.zeros([np.ceil(2*L/D).astype(int)+2,np.ceil(2*L/D).astype(int)+2])
rho = Grid(init, -L, -L, D)

#time stepping variables for animation; dt chosen through testing
dt = 0.1 #Myr
T = 200.0 #Myr, years for sun to orbit galaxy
steps = int(T/dt)

#make list of objects for plotting
images = []
#plotting setup
fig = plt.figure()
ax = plt.axes(xlim=(-L, L), ylim=(-L, L))

#assign density to grid points by inserting each array
#element from the generate galaxy function
for body in bodies:
	rho.insertBody(body)
#evaluate force on grid in particle system
rho.evalForce()
#apply force to each particle in leapfrog step
for body in bodies:
		force = rho.forceOn(body)
		body.resetForce(force[0], force[1])
		#evolve each body in time half step
		body.leapFrog(dt)
        
#evolve particle system in time
for i in range(steps):
	#counter
	print "Time step "+str(i+1)+"/"+str(steps)
        #reset density grid
        rho.resetGrid()
        #assign density to grid points by inserting each array
        #element from the generate galaxy function
        for body in bodies:
	        rho.insertBody(body)
	#evaluate force on grid in particle system
	rho.evalForce()
	#apply force to each particle
	#plt.cla()
	for body in bodies:
		force = rho.forceOn(body)
		body.resetForce(force[0], force[1])
		#evolve each body in time
		body.update(dt)
		
	#append to list of objects for plotting
	position = np.array([body.r for body in bodies]).T
	scatter, = ax.plot(position[0], position[1], 'k.')
	images.append((scatter,))

#animate plot, interval gives frames per milisecond
anim = animation.ArtistAnimation(fig, images, interval=100, blit=True)
#anim.save('ParticleMesh-Nbody'+str(N)+'.mp4
plt.show()
