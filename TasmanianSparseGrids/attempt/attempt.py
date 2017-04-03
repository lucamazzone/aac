from __future__ import division
from __future__ import print_function
import sys
sys.path.append("../InterfacePython")
sys.path.append("../../mpi4py-examples")
import TasmanianSG
import math
import struct
import testfun
import os
import numpy as np
from mpi4py import MPI
from Queue import Queue
from parutils import pprint


###################################################################
grid = TasmanianSG.TasmanianSparseGrid()

iDim = 2
iOut = 1
iDepth = 2
fTol = 1.E-5

grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, -1, "localp")
grid.setDomainTransform(np.array([[0.5, 0.8], [0.55, 0.85]]))
Points = grid.getPoints()

print(Points)


comm = MPI.COMM_WORLD

pros = sys.argv[0]

print(pros)

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
		





