#################################################################
# Name:     gradient.py                                         #
# Author:   Chris (Yuan Qi) Ni                                  #
# Date:     Sept 29, 2016                                       #
# Function: Program evaluates gradient map on images.           #
#################################################################

#essential imports
import numpy as np

#function: evaluate gradient on an image
def Grad(image, h):
    #################################################
    # image: 2D array of x, y indexed pixels        #
    # h : float pixel physical size		    #
    #################################################
    #2D arrays of partial derivatives at each point
    ddx = np.zeros(image.shape,float)
    ddy = np.zeros(image.shape,float)
    #for every row
    for i in range(image.shape[0]):
        #read every element of row
        for j in range(image.shape[1]):
            #take partial derivatives
            ddx[i][j], ddy[i][j] = Diff(image,h,i,j)
    #################################################
    # ddx: 2D array of partial x derivatives        #
    # ddy: 2D array of partial y derivatives	    #
    #################################################
    return ddx, ddy

#function: evaluate partial derivative on a mesh
def Diff(image, h, x, y):
    #################################################
    # image: 2D array of x, y indexed pixels        #
    # h : float pixel physical size		    #
    # x : int x index				    #
    # y : int y index				    #
    #################################################
    #check if pixel in on an x edge or in the center
    if x != image.shape[0]-1 and x != 0:
        #take central difference
        ddx = (image[x+1][y]-image[x-1][y])/(2*h)
    elif x == image.shape[0]-1:
        #take backward difference
        ddx = (image[x][y]-image[x-1][y])/h
    elif x == 0:
        #take forward difference
        ddx = (image[x+1][y]-image[x][y])/h
    #check if pixel in on an y edge or in the center
    if y != image.shape[1]-1 and y != 0:
        #take central difference
        ddy = (image[x][y+1]-image[x][y-1])/(2*h)
    elif y == image.shape[1]-1:
        #take backward difference
        ddy = (image[x][y]-image[x][y-1])/h
    elif y == 0:
        #take forward difference
        ddy = (image[x][y+1]-image[x][y])/h
    #################################################
    # ddx: float partial x derivative at x,y        #
    # ddy: float partial y derivative at x,y 	    #
    #################################################
    return ddx, ddy

#function: main
if __name__ == '__main__':
    #################################################
    #Sample Gradient Calculation                    #
    #################################################
    #Generate a topographical map using NASA Space  #
    #Shuttle Radar Topography Mission (SRTM) data.  #
    #################################################
    #essential imports
    import struct #for reading binary files
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    #file containing toronto topograph
    filename = "N43W080.hgt"
    f=open(filename,'rb')
    #load image from binary file (height image)
    w = np.zeros((1201,1201),float)
    #for every row
    for i in range(1201):
        #read every element of row
        for j in range(1201):
            buf = f.read(2) #read two bytes
            w[i][j] = struct.unpack('>h', buf)[0]

    #mesh distance [m]
    h = 420.0
    #Calculate partial derivatives at each point (height gradient)
    dwdx, dwdy = Grad(w, h)
    #calculate illumination intensity from height gradient
    phi = np.pi
    I = -(np.cos(phi)*dwdx + np.sin(phi)*dwdy)/np.sqrt(np.square(dwdx) + np.square(dwdy) + 1)

    #plot of topography
    surf = plt.imshow(w,vmin=70,vmax=500,cmap='gray',extent=[0,1201,0,1201])
    plt.title('Height above sea level topograph for gta')
    plt.ylabel('y [420m]')
    plt.xlabel('x [420m]')
    plt.colorbar()
    plt.show()
    #plot of intensity
    surf = plt.imshow(I,vmin=-0.01,vmax=0.01,cmap='gray',extent=[0,1201,0,1201])
    plt.title('Illumination topograph for gta (Gradient)')
    plt.ylabel('y [420m]')
    plt.xlabel('x [420m]')
    plt.colorbar()
    plt.show()