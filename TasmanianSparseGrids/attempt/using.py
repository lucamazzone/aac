from __future__ import division
from __future__ import print_function
import sys
sys.path.append("../InterfacePython")
sys.path.append("../../mpi4py-examples")
import TasmanianSG
import math
import struct
import os
from mpi4py import MPI
from queue import Queue
from operator import itemgetter
from parutils import pprint
import numpy as np
import Aggregator
import mapping


#############################################################################

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
	    
	    
WORKTAG = 1
DIETAG = 0

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
num_procs = comm.Get_size()
state = 2
loops = 7
iDim = 3
iOut = 3
iDepth = 4
fTol = 5.E-3

if my_rank == 0:
	status = MPI.Status()
	grid1 = TasmanianSG.TasmanianSparseGrid()
	grid1.makeLocalPolynomialGrid(iDim, iOut, iDepth,-1, "localp") 
	grid1.setDomainTransform(np.array([[0.7,0.81],[0.7,0.77],[0.3,0.4]]))
	Points = grid1.getPoints()
	n = len(Points)
	order = np.linspace(0,n-1,n)
	aggregate_state  = np.ones(n)*state
	Pointss = np.c_[order,Points,aggregate_state]  # qui ci vuole anche mzero , aggregate_state
	print(Pointss)
	L = []
	for i in range(n):
	    L.append([Pointss[i][0:]])
	
	wq = Work(L)
	resultz = []
	# seed slaves
	for rank in range(1,num_procs):
	    work = wq.get_next()
	    comm.send(work,dest=rank,tag=WORKTAG)
	# loop over getting new work requests until there is no more work to be done
	while True: 
	    work = wq.get_next()
	    if not work: break
	    # receive results from slave
	    result = comm.recv(source=MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status= status)
	    resultz.append(result)
	    # send slave a new work unit
	    comm.send(work, dest = status.Get_source(),tag = WORKTAG)
	# no more work to be done receive work from slaves
	for rank in range(1,num_procs):
	    result = comm.recv(source=MPI.ANY_SOURCE,tag= MPI.ANY_TAG,status=status)
	    resultz.append(result)
	# tell slaves to exit by sending empty message and DIETAG
	#for rank in range(1,num_procs):
	#    comm.send(0,dest=rank,tag=DIETAG)
	    
	results = np.vstack(resultz)
	print("after collecting")
	ff = np.zeros((n,iOut))
	for k in range(n):
	    ff[int(results[k][0])][0:] = results[k][1:]
	
	print(ff)
	
	############################################################# from here on we are actually solving
	for ciao in range(1,loops+1):
	    grid1.loadNeededPoints(ff)
	    grid1.setSurplusRefinement(fTol,-1,"fds")
	    Points = grid1.getNeededPoints()
	    aRes = grid1.evaluateBatch(Points)
	    n = len(Points)
	    order = np.linspace(0,n-1,n)
	    aggregate_state = np.ones(n)*state
	    Pointss = np.c_[order,Points,aggregate_state,aRes]
	    print(Pointss)
	    M = []
	    for i in range(n):
		    M.append([Pointss[i ][0: ]])
	    ws = Work(M)
	    resultz = []
	    for rank in range(1,num_procs):
		    work = ws.get_next()
		    comm.send(work,dest=rank, tag=WORKTAG)
		
	    while True:
		    work = ws.get_next()
		    if not work: break
		    result = comm.recv(source=MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status = status)
		    resultz.append(result)
		    comm.send(work, dest = status.Get_source(), tag = WORKTAG)
		
	    for rank in range(1,num_procs):
		    result = comm.recv(source=MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status=status)
		    resultz.append(result)
		
	    if ciao==loops:
		    for rank in range(1,num_procs):
			    comm.send(0,dest=rank,tag=DIETAG)
		    np.savetxt("Points.txt",Points)
		    np.savetxt("Predict.txt",aRes)
		    np.savetxt("Vals.txt",ff)
		
	    results = np.vstack(resultz)
	    print("after collecting, second loop")
	    ff = np.zeros((n,iOut))
	    for k in range(n):
		    ff[int(results[k][0])][0:] = results[k][1:]
		
	    #print(ff)
	
	    approx_error = np.absolute(np.subtract(aRes,ff))
	    print("mean error",np.mean(approx_error)) 
	    ff = np.add(0.75*aRes,0.25*ff)
	    
else:
	status = MPI.Status()
	while True:
	    # receive message
	    work = comm.recv(source=0, tag = MPI.ANY_TAG,status=status)
	    # check tag of received message
	    if status.Get_tag() == DIETAG: break
	    # do the work
	    resultz = np.array(work)
	    cosa = resultz[0,:].shape
	    #print(cosa[0])
	    if cosa[0] < 6: 
		    resultpp = resultz[np.ix_([0],[1,2,3])]
		    state_agg =  state        #resultz[np.ix_([0],[4])]
		    resultp = Aggregator.mapping_inverse(resultpp,state_agg)  #np.ones(3)
	    else:
		    resultpp = resultz[np.ix_([0],[1,2,3])]
		    state_agg =  state              #resultz[np.ix_([0],[4])]
		    pred = resultz[np.ix_([0],[5,6,7])]
		    (resultp,actives) = Aggregator.mapping(resultpp,state_agg,pred)
		    print(actives)
	    resulto = np.array(resultz[0][0])
	    result = np.c_[resulto,[resultp]]
	    print(result)
	    # send it back
	    comm.send(result,dest=0,tag=0)
	    
	    



#(ipsilon) = Aggregator.mapping_inverse([1.0,0.723,0.7405,0.235],1)
#print(ipsilon)

