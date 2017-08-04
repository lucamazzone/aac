import sys
#sys.path.append("../../../../users/mazzonel/aac/TasmanianSparseGrids/TSG")
sys.path.append("../InterfacePython")
#sys.path.append("../mazzonel/attempt")
import os
#script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
#rel_path = "attempt"
#abs_file_path = os.path.join(script_dir, rel_path)

#print(sys.path)


#f = open('chainsmat.txt', 'r')
#x = f.readlines()
import numpy as np
from numpy import loadtxt
import TasmanianSG
lines = loadtxt("chainsmat.txt")
#, comments="#", delimiter=",", unpack=False)

grid = TasmanianSG.TasmanianSparseGrid()
iDim = 2
iOut = 1
iDepth = 2
fTol = 1.E-5
grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, -1, "localp")
#print(lines.shape)


#print(lines[0][0]==1)
rhomat = np.array(lines)

#print(rhomat[:][11])

(a,b) = rhomat.shape
print("samples",b)
print("periods", a)
#for i in range(70):
#	if lines[1][i]==2:
#		print("ciao")

first = np.zeros(a)
second = np.zeros(a)
for i in range(a):
	first[i] = rhomat[i][0]
	second[i] = rhomat[i][1]

print(first)

order = np.linspace(0,a-1,a)
matmat = np.c_[order,first,second]


(c,d) = matmat.shape
#print(c,d)

#print(rhomat[0:][76])
#print(75.E-4)
#print(2*75E-4)

#cein = np.ones((10,5))
#cein[5][0:] = 2

#print(cein)

print(rhomat)
