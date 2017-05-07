#f = open('rhomatrix.txt', 'r')
#x = f.readlines()
import numpy as np
from numpy import loadtxt
lines = loadtxt("rhomatrix.txt")
#, comments="#", delimiter=",", unpack=False)



#print(lines.shape)

rhomat = np.array(lines)

print(rhomat)
print(rhomat[0][1])