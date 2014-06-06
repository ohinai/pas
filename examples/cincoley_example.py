


import os
pas_absolute_path = os.path.abspath("../")
import sys 
sys.path.append(pas_absolute_path)

import pas.cincoley_meng as cl
import numpy as np

case1 = cl.CincoleyMeng()

(leakage, pressure) = case1.pressure_leakage_at_t(1.)

for entry in leakage:
    print entry[0], entry[1]

