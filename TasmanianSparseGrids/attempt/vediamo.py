import numpy as np
import Aggregator
import time
start_time = time.time()
from collections import defaultdict


#Pointss = np.array([0.71,0.71,0.71])
#Altro = 1
#vediamo  = Aggregator.mapping_inverse(Pointss,Altro)
#print("--- %s seconds ---" % (time.time() - start_time))

d = {}
for x in range(1, 13):
    d['double_%02d' % x] = 'double_%02d.txt' %x

for key, value in sorted(d.items()):
    print key, value

for n in range(1,3):
	ciao = 'level_%01d.txt'  %n
	bene = np.ones(n)
	print(ciao)
	#np.savetxt(ciao,bene)
	stateofart = n*np.ones(1)
	#np.savetxt('stateofart.txt',stateofart)
	if n < 3:
		print("molto bene")


np.savetxt("stateofart.txt",0*np.ones(1))
np.savetxt("stateofart_2.txt",0*np.ones(1))
