from __future__ import division 
from __future__ import print_function 
import sys 
sys.path.append("../InterfacePython")
#sys.path.append("../../mpi4py-examples")
import TasmanianSG
import math
import struct
import os
from mpi4py import MPI
from queue import Queue
from operator import itemgetter
#from parutils import pprint
from numpy import loadtxt
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
NEWTAG = 0

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
num_procs = comm.Get_size()
loops = 9
iDim = 3
iOut = 3
iDepth = 4
fTol = 1.E-2
lowl = 0.5
highl = 0.7
lowy = 0.53
highy = 0.72
lowm = 0.05
highm = 0.4

if my_rank == 0:
	state  = 1
	status = MPI.Status() 
	grid1 = TasmanianSG.TasmanianSparseGrid()
	grid1.makeLocalPolynomialGrid(iDim, iOut, iDepth,-1, "localp")
	grid1.setDomainTransform(np.array([[lowl,highl],[lowy,highy],[lowm,highm]]))
	Points = grid1.getPoints()
	n = len(Points)
	order = np.linspace(0,n-1,n)
	aggregate_state = np.ones(n)*state
	Pointss = np.c_[order,Points,aggregate_state]
	print(Pointss)
	if NEWTAG == 1:
		np.savetxt("stateofart.txt",0*np.ones(1))
		Q = []
		for i in range(n):
	    		Q.append([Pointss[i][0:]])
	    
		wq = Work(Q)
		resultz = []
		for rank in range(1,num_procs):
	    		work = wq.get_next()
	    		comm.send(work,dest=rank,tag=WORKTAG)
		# loop getting and sending until done
		while True:
	    		work = wq.get_next()
	    		if not work: break
	    		# receive result from slave
	    		result = comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status = status)
	    		resultz.append(result)
	    		# send slave a new work unit
	    		comm.send(work,dest=status.Get_source(),tag = WORKTAG)
		# do things
		for rank in range(1,num_procs):
	    		result = comm.recv(source=MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status = status)
	    		resultz.append(result)
		# do other things
		results = np.vstack(resultz)
		ff = np.zeros((n,iOut))
		for k in range(n):
	    		ff[int(results[k][0])][0:] = results[k][1:]
		np.savetxt("ff_1.txt",ff)
		np.savetxt("stateofart.txt",0*np.ones(1))
	else:
		ff = loadtxt("ff_1.txt")
	##################################################################################### solving
	for ciao in range(1,loops+1):
	    	grid1.loadNeededPoints(ff)
	    	grid1.setSurplusRefinement(fTol,-1,"fds")
	    	fTol = fTol*1.0
	   	Points = grid1.getNeededPoints()
	    	aRes = grid1.evaluateBatch(Points)
	    	n = len(Points)
	    	order = np.linspace(0,n-1,n)
	    	aggregate_state = np.ones(n)*state
	    	Pointss = np.c_[order,Points,aggregate_state,aRes]
	    	#print(Points)
	    	print("number of points, sigma one",n)
	    	stateofart = loadtxt("stateofart.txt")
	    	if ciao>stateofart:
	    		M = []
	    		for i in range(n):
				M.append([Pointss[i][0:]])
			#	    
	    		ws = Work(M)
	    		resultz = []
	    		for rank in range(1,num_procs):
				work=ws.get_next()
		    		comm.send(work,dest=rank,tag=WORKTAG)
		    
	    		while True:
		    		work = ws.get_next()
		    		if not work: break
		    		result = comm.recv(source = MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status=status)
		    		resultz.append(result)
		    		comm.send(work, dest = status.Get_source(), tag = WORKTAG)
		    
	    		for rank in range(1,num_procs):
		    		result = comm.recv(source = MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status=status)
		    		resultz.append(result)
		    
	    		results = np.vstack(resultz)
	    		ff = np.zeros((n,iOut))
	    		for k in range(n):
		    		ff[int(results[k][0])][0:] = results[k][1:]
		    
	    		approx_error = np.absolute(np.subtract(aRes,ff))
	    		print("mean error, sigma low", np.mean(approx_error))
	    		#np.savetxt("sizeofpoints.txt",n)
	    		if ciao==loops:
				np.savetxt("Points.txt",Points)
				np.savetxt("Predict.txt",aRes)
				np.savetxt("Vals.txt",ff)
		
	    		ff = np.add(0.8*aRes,0.2*ff)
	    		titolo = 'level_%01d.txt' %ciao
			griglia = 'grid_%01d.txt' %ciao
	    		np.savetxt(titolo,ff)
			np.savetxt(griglia,Points)
			np.savetxt("stateofart.txt",ciao*np.ones(1))
		else:
			titolo = 'level_%01d.txt' %ciao
			ff = loadtxt(titolo)
	
	grid1.loadNeededPoints(ff)
	#np.savetxt("Points.txt",Points)
	#np.savetxt("Predict.txt",aRes)
	#np.savetxt("Vals.txt",ff)
################################################################################################
	state  = 2
	status = MPI.Status()
	grid2 = TasmanianSG.TasmanianSparseGrid()
	grid2.makeLocalPolynomialGrid(iDim, iOut, iDepth,-1, "localp") 
	grid2.setDomainTransform(np.array([[lowl,highl],[lowy,highy],[lowm,highm]]))
	Points = grid2.getPoints()
	n = len(Points)
	fTol = 2.E-2
	order = np.linspace(0,n-1,n)
	aggregate_state  = np.ones(n)*state
	Pointss = np.c_[order,Points,aggregate_state]  # qui ci vuole anche mzero , aggregate_state
	print(Pointss)
	if NEWTAG==1:
		np.savetxt("stateofart_2.txt",0*np.ones(1))
		L = []
		for i in range(n):
	    		L.append([Pointss[i][0:]])
	
		wq = Work(L)
		resultz = []
		# seed slaves
		for rank in range(1,num_procs):
	    		work = wq.get_next()
	    		comm.send(work,dest=rank,tag=WORKTAG)
		## loop over getting new work requests until there is no more work to be done
		while True: 
	    		work = wq.get_next()
	    		if not work: break
	    		# receive results from slave
	    		result = comm.recv(source=MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status= status)
	    		resultz.append(result)
	    		# send slave a new work unit
	    		comm.send(work,dest=status.Get_source(),tag = WORKTAG)
		# no more work to be done receive work from slaves
		for rank in range(1,num_procs):
	    		result = comm.recv(source=MPI.ANY_SOURCE,tag= MPI.ANY_TAG,status=status)
	    		resultz.append(result)
		## tell slaves to exit by sending empty message and DIETAG
		##for rank in range(1,num_procs):
		##    comm.send(0,dest=rank,tag=DIETAG)
	    
		results = np.vstack(resultz)
		ff = np.zeros((n,iOut))
		for k in range(n):
	    		ff[int(results[k][0])][0:] = results[k][1:]
	
		#print(ff)
		np.savetxt("ff_2.txt",ff)
		np.savetxt("stateofart_2.txt",0*np.ones(1))
	else:
		ff = loadtxt("ff_2.txt")
	############################################################# from here on we are actually solving
	for ciao in range(1,loops+4):
		grid2.loadNeededPoints(ff)
	    	grid2.setSurplusRefinement(fTol,-1,"fds")
	    	fTol = fTol*1.0
	    	Points = grid2.getNeededPoints()
	    	aRes = grid2.evaluateBatch(Points)
	    	n = len(Points)
	    	order = np.linspace(0,n-1,n)
	    	aggregate_state = np.ones(n)*state
	    	Pointss = np.c_[order,Points,aggregate_state,aRes]
	   	#print(Pointss)
	    	print("size of points, sigma two",n)
		stateofart_2 = loadtxt("stateofart_2.txt")
		if ciao > stateofart_2:
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
		    
	    		results = np.vstack(resultz)
	    		ff = np.zeros((n,iOut))
	    		for k in range(n):
		    		ff[int(results[k][0])][0:] = results[k][1:]
		
	    		if ciao==loops:
		    		np.savetxt("Points_2.txt",Points)
		    		np.savetxt("Predict_2.txt",aRes)
		    		np.savetxt("Vals_2.txt",ff)
	
		
	    		#print(ff)
	    		approx_error = np.absolute(np.subtract(aRes,ff))
	    		print("mean error, sigma two",np.mean(approx_error)) 
	    		ff = np.add(0.8*aRes,0.2*ff)
			titolo = 'leveltwo_%01d.txt' %ciao
			np.savetxt(titolo,ff)
			np.savetxt("stateofart_2.txt",ciao*np.ones(1))

		else:
			titolo = 'leveltwo_%01d.txt' %ciao
			ff = loadtxt(titolo)
	    
	grid2.loadNeededPoints(ff)
	##########################################################################
	############################SIMULATION####################################
	##########################################################################
	lines = loadtxt("chainsmat.txt")
	s = loadtxt("s.txt")
	chain = np.array(lines)
	#chain = np.ones((35,100))
	#chain[1][0:49] = 2
	(periods,samples) = chain.shape
	print(chain.shape)
	
	hours = 0.575*np.ones(samples)
	output = 0.61*np.ones(samples)
	meas = 0.28*np.ones(samples)
	
	hh = np.zeros((samples,periods))
	yy = np.zeros((samples,periods))
	mz = np.zeros((samples,periods))
	cc = np.zeros((samples,periods))
	hh_f = np.zeros((samples,periods))
	yy_f = np.zeros((samples,periods))
	cc_f = np.zeros((samples,periods))
	
	
	sampleorder = np.linspace(0,samples-1,samples)
	forecast = np.zeros((samples,iOut))
	Points = np.c_[hours,output,meas]
	aRes_1 = grid1.evaluateBatch(Points)
	aRes_2 = grid2.evaluateBatch(Points)
	for i in range(samples):
		if chain[0][i] ==1:
			forecast[i][0:] = aRes_1[i][0:]
		else:
			forecast[i][0:] = aRes_2[i][0:]
		
	timer = np.ones(samples)*0.0
	Sim_Points = np.c_[sampleorder,hours,output,meas,chain[0][0:],forecast,timer]
	print("siam points")
	print(Sim_Points[0][0:])
	
	for t in range(periods):
		S = []
		for j in range(samples):
		    S.append([Sim_Points[j][0:]])
	    
		ws = Work(S)
		resultz = []
		distmatrices = []
		# seed slaves again!
		for rank in range(1,num_procs):
			work = ws.get_next()
			#print(work)
			comm.send(work,dest=rank,tag=WORKTAG)
	    
		while True:
			work = ws.get_next()
			if not work: break
			result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG,status=status)
			results = result[0]
			distmatrix = result[1]
			distribution = np.dot(distmatrix,s)
			resultz.append(results)
			distmatrices.append(distribution)
			comm.send(work,dest=status.Get_source(),tag=WORKTAG)
	    
		for rank in range(1,num_procs):
			result = comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status=status)
			results = result[0]
			distmatrix = result[1]
			distribution = np.dot(distmatrix,s)
			distmatrices.append(distribution)
			resultz.append(results)
	    
		resultss = np.vstack(resultz)
		bigdistmatrix = np.vstack(distmatrices)
		print(bigdistmatrix.shape)
		fff = np.zeros((samples,iOut+1))
		big_distmatrix = np.zeros((samples,150))
		for k in range(samples):
			fff[int(resultss[k][0])][0:] = resultss[k][1:]
			big_distmatrix[int(resultss[k][0])][:] = bigdistmatrix[k][0:]
	
		#print(fff)
		for i in range(samples):
			if fff[i][1]>.9:
				hours[i] = forecast[i][0]
				output[i] = forecast[i][1]
				meas[i] = meas[i]
				hh[i][t] = forecast[i][0]
				yy[i][t] = forecast[i][1]
				mz[i][t] = meas[i]
				cc[i][t] = forecast[i][2]
			else: 
				hours[i] = fff[i][1]
				output[i] = fff[i][2]
				meas[i] = fff[i][0]
				hh[i][t] = fff[i][1]
				yy[i][t] = fff[i][2]
				mz[i][t] = fff[i][0]
				cc[i][t] = fff[i][3]
		
		timer = np.ones(samples)*(t+1)
		Points = np.c_[hours,output,meas]
	
		aRes_1 = grid1.evaluateBatch(Points)
		aRes_2 = grid2.evaluateBatch(Points)
		for i in range(samples):
			if chain[t][i]==1:
				forecast[i][0:] = aRes_1[i][0:]
			else:
				forecast[i][0:] = aRes_2[i][0:]
	
		for j in range(samples):
			hh_f[j][t] = forecast[j][0]
			yy_f[j][t] = forecast[j][1]
			cc_f[j][t] = forecast[j][2]
	
		Sim_Points = np.c_[sampleorder,hours,output,meas,chain[t][0:],forecast,timer]
		print(Sim_Points)
		np.savetxt("Hours.txt",hh)
		np.savetxt("GDP.txt",yy)
		np.savetxt("firms.txt",mz)
		np.savetxt("consumption.txt",cc)
		np.savetxt("Hours_for.txt",hh_f)
		np.savetxt("GDP_for.txt",yy_f)
		np.savetxt("cc_f.txt",cc_f)
		iteration = 'labdistrib_%01d.txt' %t
		np.savetxt(iteration,big_distmatrix)
	
	for rank in range(1,num_procs):
		comm.send(0,dest=rank,tag=DIETAG)
	###############################################################################
else:
	status = MPI.Status()
	while True:
	    # receive message
	    work = comm.recv(source=0, tag = MPI.ANY_TAG,status=status)
	    # check tag of received message
	    if status.Get_tag() == DIETAG: break
	    # do the work
	    resultz = np.array(work)
	    cosa =  max(resultz.shape)             #resultz[0,:].shape
	    #print('cosa', cosa)
	    if cosa==5:    # cosa[0] 
		    resultpp = resultz[np.ix_([0],[1,2,3])]
		    state_agg =  resultz[np.ix_([0],[4])]
		    resultp = Aggregator.mapping_inverse(resultpp,state_agg)  #np.ones(3)
	    elif cosa==8:    #cosa[0]
		    resultpp = resultz[np.ix_([0],[1,2,3])]
		    state_agg =  resultz[np.ix_([0],[4])]
		    pred = resultz[np.ix_([0],[5,6,7])]
		    (resultp,actives,momentsmat,labdist,polprime) = Aggregator.mapping(resultpp,state_agg,pred)
		    #print(actives)
	    else:
		    resultpp = resultz[np.ix_([0],[1,2,3])]
		    state_agg = resultz[np.ix_([0],[4])]
		    pred = resultz[np.ix_([0],[5,6,7])]
		    (resultp,actives,momentsmat,labdist,polprime) = Aggregator.mapping(resultpp,state_agg,pred)
		    mm  = 1.0 - actives
		    #resultp = np.c_[mm,[valls]]
		    
	    resulto = np.array(resultz[0][0])
	    if cosa==9:   # cosa[0]
		    result = []
		    results = np.c_[resulto,mm,[resultp]]
		    result.append(results)
		    result.append(labdist)
	    else: 
		    result = np.c_[resulto,[resultp]]
	    #result = np.c_[resulto,[resultp]]
	    #print(result)
	    # send it back
	    comm.send(result,dest=0,tag=0)








