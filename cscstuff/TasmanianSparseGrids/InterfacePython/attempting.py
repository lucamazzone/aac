import TasmanianSG
import numpy as np
import mapping


grid1 = TasmanianSG.TasmanianSparseGrid()

iDim = 2
iOut = 1
iDepth = 2
fTol = 1.E-5

grid1.makeLocalPolynomialGrid(iDim,iOut,iDepth,-1,"localp")

print("ciao")

points = ([1.0,2.0])

resultp = mapping.compute(points)

print(resultp)

