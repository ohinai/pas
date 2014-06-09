import numpy as np
from scipy import sparse, diag, integrate, special 
import math as math
import scipy.sparse.linalg.dsolve as dsolve
import os
from multiprocessing import Pool
import itertools

import analytical

class Thambynayagam11_2():
    
    def __init__(self):
        self.param = {
            ## Cube height in x-direction (m).
            "a":1.,
            ## Cube height in y-direction (m).
            "b":1.,
            ## Cube height in z-direction (m).
            "d":1.,
            ## Porosity.
            "porosity":1.,
            ## Compressibility (Pa^-1).
            "compressibility":1.e-8,
            ## Viscosity (Pa s).
            "viscosity":1.e-3,
            ## Permeability in the x-direction (m^2).
            "k_x":1.,
            ## Permeability in the y-direction (m^2).
            "k_y":1.,
            ## Permeability in the z-direction (m^2).
            "k_z":1.,
            ## Source term x-coordinate (m).
            "point_source_x":.5,
            ## Source term y-coordinate (m).
            "point_source_y":.5,
            ## Source term z-coordinate (m).
            "point_source_z":.5,
            }

        self.n_max = 100
        self.m_max = 100
        self.w_max = 100

    def point_source_rate(self, time):
        """ m^3/s
        """
        return 0.

    def dirichlet_0yz(self, x, y, z, t):
        return 0.

    def dirichlet_ayz(self, x, y, z, t):
        return 0.

    def dirichlet_x0z(self, x, y, z, t):
        return 0.

    def dirichlet_xbz(self, x, y, z, t):
        return 0.

    def dirichlet_xy0(self, x, y, z, t):
        return 0.

    def dirichlet_xyd(self, x, y, z, t):
        return 0.

    def point_solution(self, x, y , z, time):
        """ Returns the pressure analytical solution
        p(x, y, z, t).
        """
        n_x = self.param["k_x"]
        n_x /= self.param["porosity"]
        n_x /= self.param["compressibility"]
        n_x /= self.param["viscosity"]

        n_y = self.param["k_y"]
        n_y /= self.param["porosity"]
        n_y /= self.param["compressibility"]
        n_y /= self.param["viscosity"]

        n_z = self.param["k_z"]
        n_z /= self.param["porosity"]
        n_z /= self.param["compressibility"]
        n_z /= self.param["viscosity"]

        pressure = analytical.pointsource_all_d(x, y, z, time, 
                                                self.param["point_source_x"],
                                                self.param["point_source_y"],
                                                self.param["point_source_z"],
                                                n_x, n_y, n_z, 
                                                self.param["a"], 
                                                self.param["b"], 
                                                self.param["d"], 
                                                self.n_max, 
                                                self.m_max, 
                                                self.w_max)

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
                                            self.n_max, self.m_max, w_max)
    
    return pressure
