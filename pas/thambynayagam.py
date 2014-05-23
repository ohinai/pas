import numpy as np
from scipy import sparse, diag, integrate, special 
import math as math
import scipy.sparse.linalg.dsolve as dsolve
import os
from multiprocessing import Pool
import itertools

import analytical


def analytical_solution_d(centroid, t_in,
                          x_well, y_well, z_well, 
                          n_x, n_y, n_z, 
                          a, b, d, 
                          n_max, m_max, w_max):
    """ Computes the analytical solution for 
    all Dirichlet boundary conditions. 
    """
    (centroid, t_in,
    x_well, y_well, z_well, 
    n_x, n_y, n_z, 
    a, b, d, 
    n_max, m_max, w_max)
                          
    x_in = centroid[0]
    y_in = centroid[1]
    z_in = centroid[2]

    analytical_input = [x_in, y_in, z_in, t_in, 
                        x_well, y_well, z_well, 
                        n_x, n_y, n_z, 
                        a, b, d, 
                        n_max, m_max, w_max]

    pressure = analytical.pointsource_all_d(x_in, y_in, z_in, t_in, 
                                            x_well, y_well, z_well, 
                                            n_x, n_y, n_z, 
                                            a, b, d, 
                                            n_max, m_max, w_max)
    
    return pressure

def analytical_solution_n(centroid, t_in,
                          x_well, y_well, z_well, 
                          n_x, n_y, n_z, 
                          a, b, d, 
                          n_max, m_max, w_max):
    """ Computes the analytical solution for 
    all Neumann boundary conditions. 
    """
    x_in = centroid[0]
    y_in = centroid[1]
    z_in = centroid[2]

    analytical_input = [x_in, y_in, z_in, t_in, 
                        x_well, y_well, z_well, 
                        n_x, n_y, n_z, 
                        a, b, d, 
                        n_max, m_max, w_max]

    pressure = analytical.pointsource_all_n(x_in, y_in, z_in, t_in, 
                                            x_well, y_well, z_well, 
                                            n_x, n_y, n_z, 
                                            a, b, d, 
                                            n_max, m_max, w_max)
    
    return pressure
