#################################################################
# Name:     PMSim.py                                            #
# Authors:  Yuan Qi Ni                                          #
#           Michael Battaglia                                   #
# Date:     November 23, 2018                                   #
# Function: Model the dynamics of an N-body system with		#
#	    gravitational potentials using a particle mesh 	#
#	    method.                                             #   
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#essential imports
from MCgalaxy import generateGalaxy
from PMGrid import Grid

#Constants: Milky Way parameters
r0 = 3 #kpc, scale length of galaxy
m0 = 50.0 #10^9 solar masses, mass of galaxy
#generate 1000 masses in 15kpc box
N = 1000 #number of bodies
L = 15.0 #kpc box radius

#create bodies data
bodies = generateGalaxy(r0, m0, N, L)

#grid resolution and initializing grid
D = L/np.sqrt(N) #kpc grid spacing, based on number of bodies in our grid
init = np.zeros([np.ceil(2*L/D).astype(int)+2,np.ceil(2*L/D).astype(int)+2])
rho = Grid(init, -L, -L, D)

#time stepping variables for animation
dt = 0.1 #Myr
T = 200.0 #Myr
steps = int(T/dt)

#make list of objects for plotting
images = []
#plotting setup
fig = plt.figure()
ax = plt.axes(xlim=(-L, L), ylim=(-L, L))

#assign density to grid points by inserting bodies
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
    #assign density to grid points by inserting bodies
    for body in bodies:
        rho.insertBody(body)
    #evaluate force on grid
    rho.evalForce()
    #apply force to each particle
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
anim.save('ParticleMesh-Nbody'+str(N)+'.mp4')
plt.show()
