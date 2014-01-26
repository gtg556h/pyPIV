''' 
pyPIV version 0.10

Python implementation of Particle Image Velocimetry
Brian Williams
brn.j.williams@gmail.com
2014.01.19

'''

############################################################
''' 

Scripts required: 
    1.  Analysis for existence of feature in current frame:
        a.  Simple range?
        b.  FFT composition?
    2.  Math validation
    3.  High speed xCorr alternative (Fourier)


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
import pdb

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

        winSize = self.s0
        for i in range(0,self.nPasses):
            print('Pass',i+1)
            print('Winsize =',winSize)
            self.PIV_scan(winSize)
            winSize = np.int(np.round(winSize/2))
            
        

    #def generateGrid(self):

############################################################

    def PIV_scan(self, winSize):

        im1 = blur(self.im1,self.sigma)
        im2 = blur(self.im2,self.sigma)

        im2_padded=np.zeros([im1.shape[0]+2*winSize,im1.shape[1]+2*winSize])
        im2_padded[winSize:im2_padded.shape[0]-winSize,winSize:im2_padded.shape[1]-winSize] = im2
        
        # Analyze n-rows, m-columns:
        m = np.int(self.im1.shape[0]/winSize)
        n = np.int(self.im1.shape[1]/winSize)

        # Initialize x, y arrays for accumulating NEW shift.  Compile with self.x, self.y, only after full scan (due to potential overlap in edges)
        self.xBuild = np.zeros(im1.shape)
        self.yBuild = np.zeros(im1.shape)


        # Main grid scan:
        for i in range(0,m):
            for j in range(0,n):

                xRange = slice(i*winSize,(i+1)*winSize)
                yRange = slice(j*winSize,(j+1)*winSize)

                self.scanWindow(xRange,yRange,im1,im2_padded)


        # Address the piece on bottom row:
        if np.mod(im1.shape[0],winSize)!=0:

            for j in range(0,n):
                xRange = slice(im1.shape[0]-winSize-1, im1.shape[0])
                yRange = slice(j*winSize,(j+1)*winSize)

                self.scanWindow(xRange,yRange,im1,im2_padded)


        # Address the piece on right column:
        if np.mod(im1.shape[1],winSize)!=0:

            for i in range(0,m):
                xRange = slice(i*winSize,(i+1)*winSize)
                yRange = slice(im1.shape[1]-winSize-1,im1.shape[1])

                self.scanWindow(xRange,yRange,im1,im2_padded)


        # Update x, y arrays with current values:
        self.x = self.x + self.xBuild
        self.y = self.y + self.yBuild


        # Debugging variables:
        #self.im2_padded = im2_padded
        #self.a = a
        #self.b = b
        #self.c = c
        #self.maxIndex = maxIndex
        #self.xShift = xShift
        #self.yShift = yShift

########################################################################

    def scanWindow(self,xRange,yRange,im1,im2_padded):

        a = im1[xRange, yRange]

            # Add check: if range[b] < min, skip cross_correlation
            # (To account for window with no particles!)

            if 1:

                xShift = self.x[np.int(np.mean([xRange.start,xRange.stop])),np.int(np.mean([yRange.start,yRange.stop]))]
                yShift = self.y[np.int(np.mean([xRange.start,xRange.stop])),np.int(np.mean([yRange.start,yRange.stop]))]

                b = im2_padded[xRange.start + xShift:xRange.end + 3*winSize + xShift, yRange.start + yShift:yRange.end + 3*winSize + yShift]

                c, maxIndex,xCurrent,yCurrent = xcorr(a,b,winSize)

                # Note: in following lines, I'm overwriting some previously computed in the central region during border scans...  Probably not an issue...
                self.xBuild[xRange,yRange] = xCurrent
                self.yBuild[xRange,yRange] = yCurrent



############################################################
############################################################
############################################################
############################################################

def xcorr(a,b,winSize):
    c = scipy.signal.correlate2d(a.astype(float),b.astype(float))
    max_in = np.unravel_index(c.argmax(),c.shape)
    xCurrent = 4*winSize/2 - 1 - max_in[0]
    yCurrent = 4*winSize/2 - 1 - max_in[1]
    print("x=",xCurrent)
    print("y=",yCurrent)

    return c, max_in, xCurrent, yCurrent
    
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
    
