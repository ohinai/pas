


import os
pas_absolute_path = os.path.abspath("../")
import sys 
sys.path.append(pas_absolute_path)

import pas.cincoley_meng as cl
import numpy as np


(leakage, pressure) = cl.cinco_ley_laplace(5., 10.)

for (x, f) in enumerate(leakage):
    print x, f


