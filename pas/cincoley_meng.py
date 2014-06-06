

import numpy as np
from scipy import sparse, diag, integrate, special 
import math as math
import scipy.sparse.linalg.dsolve as dsolve
import os
from multiprocessing import Pool
import itertools

def gaver_stehfest(time, lap_func):
    """ Performs a numerical inversion
    of the Laplace transform using the
    Gaver-Stehfest algorithm. The algorithm
    can be found in:
    "A Uninfied Framework for Numerically Inverting Laplace
    Transforms" By Joseph Abate and Ward Whitt.
    """
    def nCr(n, r):
        return math.factorial(n)/(math.factorial(r)*
                                  math.factorial(n-r))
    def a(k, n):
        summation = 0.
        for j in range((k+1)/2, min(k, n)+1):
            current_summation = float(pow(j, n+1))/float(math.factorial(n))
            current_summation *= nCr(n, j)
            current_summation *= nCr(2*j, j)
            current_summation *= nCr(j, k-j)
            summation += current_summation
        return summation*pow(-1, n+k)
    n = 7
    total_sum = a(1, n)*lap_func(1.*np.log(2.)/time)
    for k in range(2, 2*n+1):
        total_sum += a(k, n)*lap_func(k*np.log(2.)/time)
    return total_sum*np.log(2.)/time
    
class CincoleyMeng():
    """ Computes the semi-anlytical solution from 
    Cinco Ley and Meng:
    "Pressure Transient Analysis of Wells With Finite Conductivity 
    Vertical Fractures in Double Porosity Reservoirs"
    Also see:
    "Non-Darcy Flow in Wells With Finite-Conductivity 
    Vertical Fractures. 

    :ivar float c_fd: Fracture conductivity difined as \
    (fracture_k*fracture_width)/(reservoir_k*fracture_length). 
    
    t_D = \frac{\Beta res_k t}{\theta c_t \mu x_f^2}
    p_D = \frac{res_k h [p_i - p(t)]}{\alpha q \Beta \mu}
    q_D = \frac{2 q_f}{q_w}
    """
    def __init__(self):

        self.param = {
            ## Reservoir permeability (m^2).
            "res_k":1.e-3,
            ## Fracture permeability. 
            "frac_k":1.,
            ## Fracture length (m).
            "frac_length":1.,
            ## Fracture width (m). 
            "frac_width":.01,
            ## Fluid viscosity (Pa s). 
            "viscosity":1.e-3, 
            ## Reservoir porosity. 
            "res_porosity":1.,
            ## Compressibility (Pa^-1). 
            "compressibility":1.e-8}

        self.number_of_segments = 40

    def pressure_leakage_at_t(self, time):
        """ Computes fracture leakage and bottom hole pressure at non-dimensional time t_d. 
        
        :param float t_d: Non-dimensional time. 
        :return (numpy.array, float): Returns a tuple, with the first entry \
        being the leakage over the length of the fracture. The \
        number of entries corresponds to the number_of_segments used for the
        semi-analytical solution. The second entry in the tuple is the bottom 
        hole pressure at that time. 
        """
        t_d = self.param["res_k"]*time
        t_d /= self.param["res_porosity"]
        t_d /= self.param["viscosity"]
        t_d /= self.param["compressibility"]
        t_d /= self.param["frac_length"]**2

        x_d = np.array([float(i)/self.number_of_segments for i in range(self.number_of_segments)])

        c_fd = self.param["frac_k"]
        c_fd *= self.param["frac_width"]
        c_fd /= self.param["res_k"]
        c_fd /= self.param["frac_length"]
        
        delta_x = 1./float(self.number_of_segments)

        def laplace_solution(s):
            lhs = np.zeros((self.number_of_segments+1, self.number_of_segments+1))
            rhs = np.zeros(self.number_of_segments+1)
            for j in range(self.number_of_segments):
                x_dj = x_d[j]+.5*delta_x
                rhs[j] = np.pi*x_dj/(c_fd*s)
                for i in range(self.number_of_segments):
                    x_di = x_d[i]
                    if x_dj > x_di and x_dj < x_di+delta_x:
                        (quad1, error1) = integrate.quad(lambda x: special.kn(0, abs(x_dj-x)*np.sqrt(s)), 
                                                         x_di,
                                                         x_di+delta_x,
                                                         points = [x_dj], epsrel=1.e-16)
                    else:
                        (quad1, error1) = integrate.quad(lambda x: special.kn(0, abs(x_dj-x)*np.sqrt(s)), 
                                                         x_di,
                                                         x_di+delta_x, epsrel=1.e-16)

                    (quad2, error2) = integrate.quad(lambda x: special.kn(0, abs(x_dj+x)*np.sqrt(s)), 
                                                     x_di,
                                                     x_di+delta_x, epsrel=1.e-16)
                    #assert (error1 < 1.e-8)
                    #assert (error2 < 1.e-8)
                    lhs[j, i] += -.5*(quad1+quad2)
                    if i < j:
                        lhs[j, i] += (np.pi/c_fd)*((delta_x**2)/2.+
                                                   (x_dj-(i+1)*delta_x)*delta_x)

                lhs[j, j] += (np.pi/c_fd)*(delta_x**2)/8.

            for i in range(self.number_of_segments):
                lhs[i, self.number_of_segments] += 1.

            for j in range(self.number_of_segments):
                lhs[self.number_of_segments, j] += delta_x

            rhs[self.number_of_segments] = 1./s
            solution = np.linalg.solve(lhs, rhs)
            min_singular_value = np.linalg.svd(lhs)[1][-1]
            assert(min_singular_value > 1.e-8)

            return solution

        full_solution = gaver_stehfest(t_d, laplace_solution)

        return (zip(x_d*self.param["frac_length"], full_solution[:-1]), full_solution[-1])

