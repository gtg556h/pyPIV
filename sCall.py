import numpy as np
import matplotlib.pyplot as plt
import pivLib


file1 = 'synt0.png'
file2 = 'synt1.png'
show = 0
sigma = 2
s0 = 50
sRatio = 0.5
nPasses = 3

params = {'file1':file1, 'file2':file2, 'show':show, 's0':s0, 'sRatio':sRatio, 'nPasses':nPasses, 'sigma':sigma}

c = pivLib.pivEval(params)




