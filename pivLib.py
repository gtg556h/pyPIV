import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
from scipy import misc


class pivEval(object):

    def __init__(self, params):
        #print Initializing
        self.im1 = misc.imread(params['file1'])
        self.im2 = misc.imread(params['file2'])

        if params['show'] == 1:
            plt.subplot(121)
            plt.imshow(self.im1)
            plt.subplot(122)
            plt.imshow(self.im2)
            plt.show()



