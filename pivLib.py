''' 
pyPIV version 0.10

Python implementation of Particle Image Velocimetry
Brian Williams
brn.j.williams@gmail.com
2014.01.19

'''

############################################################
############################################################

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy.linalg
from scipy import misc 
import scipy.signal
import scipy.ndimage

############################################################
############################################################

class pivEval(object):

    def __init__(self, params):
        #print Initializing
        self.im1 = misc.imread(params['file1'])
        self.im2 = misc.imread(params['file2'])
        self.s0 = params['s0']
        self.nPasses = params['nPasses']
        self.sRatio = params['sRatio']
        self.sigma = params['sigma']

        self.x = np.zeros(np.shape(self.im1))
        self.y = np.zeros(np.shape(self.im1))

        if params['show'] == 1:
            plt.subplot(121)
            plt.imshow(self.im1)
            plt.subplot(122)
            plt.imshow(self.im2)
            plt.show()

        self.PIV_scan(50)
        

    #def generateGrid(self):


    def PIV_scan(self, winSize):

        im1 = blur(self.im1,self.sigma)
        im2 = blur(self.im2,self.sigma)


        im2_padded=np.zeros([im1.shape[0]+2*winSize,im1.shape[1]+2*winSize])
        im2_padded[winSize:im2_padded.shape[0]-winSize,winSize:im2_padded.shape[1]-winSize] = im2
        for i in range(0,1):#range(0,np.int(self.im1.shape[0]/winSize)):
            for j in range(0,1):#range(0,np.int(self.im1.shape[1]/winSize)):
                a = im1[i*winSize:(i+1)*winSize,j*winSize:(j+1)*winSize]
                b = im2_padded[i*winSize:(i+3)*winSize,j*winSize:(j+3)*winSize]
                c, maxIndex = xcorr(a,b)
                self.x[i*winSize:(i+1)*winSize,j*winSize:(j+1)*winSize] = self.x[i+winSize:(i+1)*winSize,j*winSize:(j+1)*winSize] # + !!!!!!!!!!!!! Compute!!!

        self.im2_padded = im2_padded
        self.a = a
        self.b = b
        self.c = c
        self.maxIndex = maxIndex



############################################################
############################################################
############################################################
############################################################

def xcorr(a,b):
    c = scipy.signal.correlate2d(a.astype(float),b.astype(float))
    max_in = np.unravel_index(c.argmax(),c.shape)
    return c, max_in
    
############################################################

def blur(z, sigma=2):
    blurred = scipy.ndimage.gaussian_filter(z, sigma=sigma)
    return blurred

############################################################

def plotXCorr(Z):
    X = np.arange(0,Z.shape[0])
    Y = np.arange(0,Z.shape[1])
    X, Y = np.meshgrid(X, Y)

    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01,1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
    
