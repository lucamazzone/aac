import TasmanianSG
import numpy as np



grid1 = TasmanianSG.TasmanianSparseGrid()

iDim = 2
iOut = 1
iDepth = 2
fTol = 1.E-5

grid1.makeLocalPolynomialGrid(iDim,iOut,iDepth,-1,"localp")

