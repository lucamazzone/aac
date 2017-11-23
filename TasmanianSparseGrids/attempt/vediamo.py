import numpy as np
import Aggregator
import Stuff
import time
start_time = time.time()
from collections import defaultdict
from numpy import loadtxt

#Pointss = np.array([0.71,0.71,0.71])
#Altro = 1
#vediamo  = Aggregator.mapping_inverse(Pointss,Altro)
#print("--- %s seconds ---" % (time.time() - start_time))

d = {}
for x in range(1, 13):
    d['double_%02d' % x] = 'double_%02d.txt' %x

#for key, value in sorted(d.items()):
#    print key, value

#for n in range(1,3):
#	ciao = 'level_%01d.txt'  %n
#	bene = np.ones(n)
#	print(ciao)
#	#np.savetxt(ciao,bene)
#	stateofart = n*np.ones(1)
#	#np.savetxt('stateofart.txt',stateofart)
#	if n < 3:
#		print("molto bene")


#np.savetxt("stateofart.txt",0*np.ones(1))
#np.savetxt("stateofart_2.txt",0*np.ones(1))


resultpp = np.array([0.6,0.65,0.35])
pred = np.array([0.58,0.55,0.63])
state_agg = 1

#(resultp,actives,momentsmat1,labdist1,polprime1) = Aggregator.mapping(resultpp,state_agg,pred)

#print(resultp)
#print(actives)
#print(momentsmat1)

#resultpp = np.array([0.6,0.65,0.35])
#pred = np.array([0.57,0.57,0.63])
#state_agg = 1

#(resultp,actives,momentsmat2,labdist2,polprime2) = Aggregator.mapping(resultpp,state_agg,pred)

#print(resultp)
#print(momentsmat2)

#prova = np.dstack((momentsmat1, momentsmat2))
#print(prova.shape)

#bumbum = np.sum(prova,axis=2)/2 #prova[:,:,0]+prova[:,:,1]

#list3=[]
#list3.append(resultp)
#list3.append(bumbum)

#print(list3)
#print(np.array(list3[0]))
#print(np.array(list3[1]))

#prova2 = np.empty_like(prova)
#prova2[:,:,0] = prova[:,:,1]
#prova2[:,:,1] = prova[:,:,0]
#print(prova2.shape)

#print(prova2[:,:,0])


list4 = []
#list4.append(momentsmat1)
#list4.append(momentsmat2)

#bohm = np.dstack(list4)

#print(bohm.shape)
#print(bohm[:,:,1])

s = loadtxt("s.txt")

print(s)

#print(np.dot(s,momentsmat1))
#print(np.dot(s,momentsmat2))


#rhomat = loadtxt("rhomatrix.txt")

#zct = 7
#kval = 0.3
#agg_state = 1
#momnum = 5

#rhomatt = rhomat[:,agg_state-1]
#rhomatt = rhomatt.reshape(5,10).T
#print(rhomatt)
#rhovec = rhomatt[zct,:]
#(weights,nodes) = Aggregator.qsimpweightsnodes(0.1,1.3,50)

#momvec = momentsmat1[zct,:]


#print(momvec)
#print(rhovec)

#F_k = rhovec[0]*(kval- momvec[0])

#for j in range(1,momnum):
#	F_k = F_k + rhovec[j]*(  (kval - momvec[0])**j - momvec[j] )

#F_k = np.exp(F_k)



#(weights,nodes) = Aggregator.qsimpweightsnodes(0.1,1.3,50)


#distribution = np.dot(labdist2,s)
#print(distribution)

chainsmat = loadtxt("chainsmat.txt")




#chainsmat2 = chainsmat

#chainsmat[25:26,0:59] = 2

#print(chainsmat.shape)

#chainsmat = np.hstack((chainsmat,chainsmat2))

print(chainsmat.shape)
print(chainsmat[25,0:59])
print(chainsmat[25,60:119])
print(chainsmat[26,0:59])
print(chainsmat[26,60:119])

chainsmat[25,60:119] = 1


np.savetxt("chainsmat.txt",chainsmat)


s1 = np.random.uniform(0.5,0.7,20)
s2 = np.random.uniform(0.52,0.72,20)
s3 = np.random.uniform(0.10,0.35,20)

bigx = np.vstack((s1,s2,s3)).T

#print(bigx)

n = len(bigx)

#print(n)

#for k in range(n):
#	print(k)
