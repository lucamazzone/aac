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


###################################################################

class Work(object):
    def __init__(self,points):
	    q = Queue()
	    for p in points:
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

if my_rank == 0:
	status = MPI.Status()
	grid = TasmanianSG.TasmanianSparseGrid()
	iDim = 2
	iOut = 1
	iDepth = 2
	fTol = 1.E-5
	grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, -1, "localp")
	grid.setDomainTransform(np.array([[0.5, 0.8], [0.55, 0.85]]))
	Points = grid.getPoints()
	n = len(Points)
	print(Points)
	order = np.linspace(0,n-1,n)
	Pointss = np.c_[order, Points]
	print(Pointss)
	L =  []         #[[Points[0][0],Points[0][1]]]
	for i in range(n):
	    L.append([Pointss[i][0:]])
	
	wq = Work(L)
	#work = wq.get_next()
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
	for rank in range(1,num_procs):
	    comm.send(0,dest=rank,tag=DIETAG)
	
	results = np.vstack(resultz)  #array
	print("after collecting")
	ff = np.zeros((n,2))
	for k in range(n):
		ff[int(results[k][0])][0:] = results[k][1:]
	print(ff)
	print(wq.get_next())
	ws = Work(L)
	print(ws.get_next())
else:
	print(my_rank)
	status = MPI.Status()
	while True:
	    # receive message from master
	    work = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
	    # check tag of received message
	    if status.Get_tag() == DIETAG: break
	    # do the work
	    resultz = np.array(work)
	    resultp =2*resultz[np.ix_([0],[1,2])]
	    resulto = np.array(resultz[0][0])
	    result = np.c_[resulto,resultp]
	    print(result)
	    # send results back
	    comm.send(result,dest=0,tag=0)




#q = Queue()

#print(q)


#pprint("-"*78)
#pprint(" Running on %d cores" % comm.size)
#pprint("."*78)

#rank = MPI.COMM_WORLD.Get_rank()
# print("[%d]" % rank)

#if comm.rank == 0:
#    grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, -1, "localp") 
#    Points = grid.getPoints()
#    aVals = np.empty([Points.shape[0],2])
#else:
#    Points = np.empty([13,2])
#    aVals = np.empty([Points.shape[0],2])


#rec_points = np.empty([1,2])
#rec_vals = np.empty([1,2])

# scatter
#comm.Scatter([Points, MPI.DOUBLE],[rec_points, MPI.DOUBLE])

#pprint("After Scatter")
#for r in range(comm.size):
#	if MPI.COMM_WORLD.Get_rank() == r:
#	    pprint("[%d] of [%s]" %(rank,comm.size) )    
#	comm.Barrier()
#print(rec_points)


#rec_vals[0,1] = 1.0
#rec_vals[0,0] = testfun.funfun(rec_points)
 
#print(rec_vals)

#comm.Allgather([rec_vals, MPI.DOUBLE] ,[aVals , MPI.DOUBLE])

#pprint("After Gather")
#for r in range(comm.size):
#	if comm.rank == int(r):
#pprint("%d" % comm.rank)
# print(aVals)
		





