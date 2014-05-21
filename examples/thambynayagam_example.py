
import pas.thambynayagam as tham
import numpy as np

k = 1.e-12
porosity = .2
viscosity = 1.e-4
compressibility = 1.e-8

n_x = k
n_x /= porosity
n_x /= viscosity
n_x /= compressibility

n_y = n_x
n_z = n_x 

a = 1
b = 1
d = 1

n_max = 100
m_max = 100
w_max = 100

pressure = tham.analytical_solution_d(np.array([.1, .1, .1]), 1.,
                                      .5, .5, .5, 
                                      n_x, n_y, n_z, 
                                      a, b, d, 
                                      n_max, m_max, w_max)

print pressure
