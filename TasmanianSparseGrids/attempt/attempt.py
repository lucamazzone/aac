from __future__ import division
from __future__ import print_function
import sys
sys.path.append("../InterfacePython")
sys.path.append("../../mpi4py-examples")
import TasmanianSG
import math
import struct
import os
import numpy as np
from mpi4py import MPI
from queue import Queue
from operator import itemgetter
from parutils import pprint
from numpy import f2py
import mapping
#import matplotlib.pyplot as plt
    


###################################################################

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
loops = 2

if my_rank == 0:
	status = MPI.Status()
	grid = TasmanianSG.TasmanianSparseGrid()
	iDim = 2
	iOut = 1
	iDepth = 2
	fTol = 1.E-5
	grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, -1, "localp")
	#grid.setDomainTransform(np.array([[0.5, 0.8], [0.55, 0.85]]))
	Points = grid.getPoints()
	n = len(Points)
	order = np.linspace(0,n-1,n)
	Pointss = np.c_[order, Points]
	print(Pointss)
	L =  []       
	for i in range(n):
	    L.append([Pointss[i][0:]])
	
	wq = Work(L)
	resultz = []
	# seed slaves and send them one unit of work each
	for rank in range(1,num_procs):
	    work = wq.get_next()
	    comm.send(work,dest=rank, tag=WORKTAG)
	
	# loop over getting new work requests until there is no more work to be done
	while True:
	    work = wq.get_next()
	    if not work: break
	    # receive results from slave
	    result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
	    resultz.append(result)
	    # send slave a new work unit
	    comm.send(work, dest = status.Get_source(), tag=WORKTAG)
	# no morework to be done, receive work from slaves
	for rank in range(1,num_procs):
	    result = comm.recv(source=MPI.ANY_SOURCE, tag =MPI.ANY_TAG, status=status)
	    resultz.append(result)
	# tell slaves to exit by sending empty message with DIETAG
#	for rank in range(1,num_procs):
#	    comm.send(0,dest=rank,tag=DIETAG)
	#cosa = result[0,:].shape
	#print("size of result", cosa[0])
	#if cosa[0] > 1:
	#	print("MOOOOLTO BENE")
	
	results = np.vstack(resultz)  #array
	print("after collecting")
	ff = np.zeros((n,iOut))
	for k in range(n):
		ff[int(results[k][0])][0:] = results[k][1:]
	print(ff)
	##############################################################
	for ciao in range(1,loops+1):
	    grid.loadNeededPoints(ff)
	    grid.setSurplusRefinement(fTol,-1,"fds")
	    Points = grid.getNeededPoints()
	    n = len(Points)
	    order = np.linspace(0,n-1,n)
	    Pointss = np.c_[order,Points]
	    print(Pointss)
	    M = []
	    for i in range(n):
		    M.append([Pointss[i][0:]])
	    ws = Work(M)
	    resultz = []
	    for rank in range(1,num_procs):
		    work = ws.get_next()
		    comm.send(work,dest=rank, tag=WORKTAG)
		
	    while True:
		    work = ws.get_next()
		    if not work: break
		    result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
		    resultz.append(result)
		    comm.send(work, dest = status.Get_source(), tag=WORKTAG)
	
	    for rank in range(1,num_procs):
		    result = comm.recv(source=MPI.ANY_SOURCE, tag =MPI.ANY_TAG, status=status)
		    resultz.append(result)
		
	    print(ciao)
	    if ciao==loops:
		    print("loops", ciao)
		    for rank in range(1,num_procs):
			    comm.send(0,dest=rank,tag=DIETAG)
	
	    results = np.vstack(resultz)
	    print("after collecting, second loop")
	    ff = np.zeros((n,iOut))
	    for k in range(n):
		    ff[int(results[k][0])][0:] = results[k][1:]
	
	    print(ff)	
	    aRes = grid.evaluateBatch(Points)
	    error =  np.absolute(np.subtract(aRes,ff))
	    print("mean abs error is", np.mean(error))
	#######################################################################
#	grid.loadNeededPoints(ff)
#	grid.setSurplusRefinement(fTol,-1,"fds")
#	Points = grid.getNeededPoints()
#	n = len(Points)
#	order = np.linspace(0,n-1,n)
#	Pointss = np.c_[order,Points]
#	print(Pointss)
#	N = []
#	for i in range(n):
#		N.append([Pointss[i][0:]])
#	wr = Work(N)
#	resultz = []
#	for rank in range(1,num_procs):
#		work = wr.get_next()
#		comm.send(work,dest=rank, tag= WORKTAG)
#		
#	while True:
#		work = wr.get_next()
#		if not work: break
#		result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
#		resultz.append(result)
#		comm.send(work,dest = status.Get_source(), tag = WORKTAG)
#		
#	for rank in range(1,num_procs):
#		result = comm.recv(source=MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status=status)
#		resultz.append(result)
#		
#	for rank in range(1,num_procs):
#		comm.send(0,dest=rank,tag=DIETAG)
#		
#	results = np.vstack(resultz)
#	print("after collecting, third loop")
#	ff = np.zeros((n,iOut))
#	for k in range(n):
#		ff[int(results[k][0])][0:] = results[k][1:]
#		
#	print(ff)
#	aRes = grid.evaluateBatch(Points)
#	error = np.absolute(np.subtract(aRes,ff))
#	print("mean abs error is", np.mean(error))
	##################################################################################
else:
	status = MPI.Status()
	while True:
	    # receive message from master
	    work = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
	    # check tag of received message
	    if status.Get_tag() == DIETAG: break
	    # do the work
	    resultz = np.array(work)
	    resultpp = resultz[np.ix_([0],[1,2])]
	    resultp = mapping.compute(resultpp)
	    print("resultp",resultp)
	    resulto = np.array(resultz[0][0])
	    result = np.c_[resulto,[resultp]]
	    print(result)
	    # send results back
	    comm.send(result,dest=0,tag=0)


#plt.scatter(Points[0][0:], Points[1][0:])
#plt.show()









