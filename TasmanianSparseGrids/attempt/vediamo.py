import numpy as np
import Aggregator
import time
start_time = time.time()

Pointss = np.array([0.71,0.71,0.71])
Altro = 1
vediamo  = Aggregator.mapping_inverse(Pointss,Altro)
print("--- %s seconds ---" % (time.time() - start_time))

