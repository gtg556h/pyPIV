import numpy as np
import matplotlib.pyplot as plt
import pivLib
import scipy.misc
#import ipdb
import pdb


a = np.zeros([50,50])
b = np.zeros([50,50])
a[25,25] = 10
a[10,10] = 10
a[10,40] = 10
a[40,10] = 10
a[40,40] = 10
b[30,25] = 10
b[15,10] = 10
b[15,40] = 10
b[45,10] = 10
b[45,40] = 10
a = pivLib.blur(a,2)
b = pivLib.blur(b,2)
scipy.misc.imsave('a.png',a)
scipy.misc.imsave('b.png',b)


file1 = 'a.png'
file2 = 'b.png'
show = 0
sigma = 2
s0 = 50
sRatio = 0.5
nPasses = 2

params = {'file1':file1, 'file2':file2, 'show':show, 's0':s0, 'sRatio':sRatio, 'nPasses':nPasses, 'sigma':sigma}

c = pivLib.pivEval(params)



