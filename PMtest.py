#################################################################
# Name:     PMtest.py                                           #
# Authors:  Yuan Qi Ni                                          #
#           Michael Battaglia                                   #
# Date:     November 23, 2018                                   #
# Function: Model the dynamics of an N-body system with		#
#           gravitational potentials using a particle mesh 	#
#	    method. Testing algorithms.          	        #   
#################################################################

#essential modules
import time
import numpy as np
import matplotlib.pyplot as plt

#essential imports
from PMGrid import Grid
from MCgalaxy import generateGalaxy

#function: main
if __name__ == '__main__':
    #################################
    # Computational Complexity Test #
    #################################
    
    #Constants: Milky Way parameters
    r0 = 3 #kiloparsecs, scale length of galaxy
    m0 = 50.0 #10^9 solar masses, mass of galaxy
    #generate 1000 masses in 15kpc box
    N = 1000 #number of bodies
    L = 15.0 #Number of grid divisions, decided on after trying several divisions. 
    
    #Create bodies data
    bodies = generateGalaxy(r0, m0, N, L)
    
    #grid resolution and initializing grid
    D = L/np.sqrt(10*N) #kpc grid spacing, based on number of bodies in our grid
    init = np.zeros([np.ceil(2*L/D).astype(int)+2,np.ceil(2*L/D).astype(int)+2])
    rho = Grid(init, -L, -L, D)
    
    #time stepping variables for animation
    dt = 0.1 #Myr
    T = 200.0 #Myr
    steps = int(T/dt)
    
    #list of N for which to test runtime
    elapsed_list = []
    N_list = range(1, 101)+range(101,1001,10)+range(1001,10001,100)
    #calculating time to evaluate mesh for each N
    for i in N_list:
        #generate bodies
        i_bodies = generateGalaxy(r0, m0, i, L)
        #start counter
        start = time.clock()
        #add bodies to grid
        for i_body in i_bodies:
            rho.insertBody(i_body)
        #evaluate force on grid in particle system
        rho.evalForce()
        #apply force to each particle
        for i_body in i_bodies:
            force = rho.forceOn(i_body)
            i_body.resetForce(force[0], force[1])
            #evolve each body in time
            i_body.update(dt)
        elapsed = (time.clock() - start)
        elapsed_list.append(elapsed)
        print "N step "+str(i)+"/10000"
    #plot runtime
    plt.plot(N_list, elapsed_list,label= 'PM Runtime')
    #NlogN time
    NlnN = N_list*np.log(N_list)
    #normalize
    NlnN = NlnN*(elapsed_list[-1]/NlnN[-1])
    #plot NlogN time
    plt.plot(N_list, NlnN, label="NlogN")
    plt.xlabel("Number of bodies")
    plt.ylabel("Time elapsed (s)")
    plt.title("Computational time for increasing N")
    plt.legend()
    plt.show()

    ############################
    # Energy Conservation Test #
    ############################
    
    #Create bodies data
    bodies = generateGalaxy(r0, m0, N, L)

    #grid resolution and initializing grid
    D = L/np.sqrt(N) #kpc grid spacing, based on number of bodies in our grid
    init = np.zeros([np.ceil(2*L/D).astype(int)+2,np.ceil(2*L/D).astype(int)+2])
    rho = Grid(init, -L, -L, D)

    #time stepping variables for animation
    dt = 0.1 #Myr
    T = 200.0 #Myr
    steps = int(T/dt)

    #assign density to grid points by inserting bodies
    for body in bodies:
        rho.insertBody(body)

        '''
        #Testing function. Plots densities on the grid
        #with weight one.
        #body.plot()
    #rho.plot(m0/N)
    #plt.xlim([-15,15])
    #plt.ylim([-15,15])
    '''
    #evaluate force on grid in particle system
    rho.evalForce()
    #apply force to each particle in leapfrog step
    for body in bodies:
        force = rho.forceOn(body)
        body.resetForce(force[0], force[1])
        #evolve each body in time half step
        body.leapFrog(dt)
        
    #track energy at each time step.
    t = np.linspace(0,steps*dt,steps)
    E = np.zeros(len(t))
    #sum kinetic
    for body in bodies:
        E[0]+=body.Kenergy(dt)
    #sum potential
    for j in range(len(bodies)):
        for k in range(j+1,len(bodies)):
            E[0]+=bodies[j].Uinteract(bodies[k])

    #evolve particle system in time
    #WARNING, slow because O(N^2) calculation
    for i in range(steps):
        #counter
        print "Time step "+str(i+1)+"/"+str(steps)
        #evaluate force on grid in particle system
        rho.evalForce()
        #apply force to each particle
        for body in bodies:
            force = rho.forceOn(body)
            body.resetForce(force[0], force[1])
            #evolve each body in time
            body.update(dt)
        #calculate energy at each time step
        for body in bodies:
            E[i]+=body.Kenergy(dt)
        for j in range(len(bodies)):
            for k in range(j+1,len(bodies)):
                E[i]+=bodies[j].Uinteract(bodies[k])

    #Energy plot
    plt.plot(t, E)
    plt.title("Energy conservation")
    plt.ylabel("Energy"+r'[$M\odot kpc^2/(10Myr)^2$]')
    plt.show()
