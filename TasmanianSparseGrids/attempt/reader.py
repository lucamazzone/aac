#f = open('rhomatrix.txt', 'r')
#x = f.readlines()
import numpy as np
from numpy import loadtxt
lines = loadtxt("chainsmat.txt")
#, comments="#", delimiter=",", unpack=False)



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



order = np.linspace(0,a-1,a)
matmat = np.c_[order,first,second]


(c,d) = matmat.shape
#print(c,d)

#print(rhomat[0:][76])
#print(75.E-4)
#print(2*75E-4)

cein = np.ones((10,5))
cein[5][0:] = 2

print(cein)

