import numpy as np
from numpy import loadtxt
from queue import Queue


class Work(object):
	def __init__(self,inputs):
		q = Queue()
		for p in inputs:
			q.put(p)
		self.work = q

	def get_next(self):
		if self.work.empty():
			return None
		return self.work.get()





iOut = 3
lines = loadtxt("chainsmat.txt")
chain = np.array(lines)
(periods,samples) = chain.shape
print(chain.shape)
print(periods)


#print(chain)


hours = 0.7*np.ones(samples)
output=0.75*np.ones(samples)
meas = 0.22*np.ones(samples)
sampleorder = np.linspace(0,samples-1,samples)
forecast = np.zeros((samples,iOut))
Points = np.c_[hours,output,meas]
aRes_1 = Points
aRes_2 = Points*2.0

print(aRes_1)
print(Points.shape)

for i in range(samples):
	if chain[0][i]==1:
		forecast[i][0:] = aRes_1[i][0:]
	else:
		forecast[i][0:] = aRes_2[i][0:]

timer = np.ones(samples)*0.0

Sim_Points = np.c_[sampleorder,Points,chain[0][0:],forecast,timer]


print(Sim_Points.shape)
print(Sim_Points)


S = []

for j in range(samples):
	S.append([Sim_Points[j][0:]])


#print(S)
ws = Work(S)
resultz = []


for rank in range(1,12):
	work = ws.get_next()
	print(work)
