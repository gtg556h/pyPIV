import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy.linalg
from scipy import misc 
import scipy.signal
import scipy.ndimage


class pivEval(object):

    def __init__(self, params):
        #print Initializing
        self.im1 = misc.imread(params['file1'])
        self.im2 = misc.imread(params['file2'])
        self.s0 = params['s0']
        self.nPasses = params['nPasses']
        self.sRatio = params['sRatio']

        if params['show'] == 1:
            plt.subplot(121)
            plt.imshow(self.im1)
            plt.subplot(122)
            plt.imshow(self.im2)
            plt.show()



def xcorr(a,b):
    c = scipy.signal.correlate2d(a.astype(float),b.astype(float))
    return c


def blur(z, sigma=3):
    blurred = scipy.ndimage.gaussian_filter(z, sigma=sigma)
    return blurred


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
    
