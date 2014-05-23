

import numpy as np
from scipy import sparse, diag, integrate, special 
import math as math
import scipy.sparse.linalg.dsolve as dsolve
import os
from multiprocessing import Pool
import itertools

def gaver_stehfest(t, P):
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
    total_sum = a(1, n)*P(1.*np.log(2.)/t)
    for k in range(2, 2*n+1):
        total_sum += a(k, n)*P(k*np.log(2.)/t)
    return total_sum*np.log(2.)/t
    
def CincoleyMeng():
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
        
        self.c_fd = None
        self.fracture_length = 1.
        self.fracture_width = .01
        self.fracture_k = 1.
        self.reservoir_k = 1.
        self.viscosity = 1000.
        
        self.number_of_segments = 20

    def pressure_leakage_at_t_d(self, t_d):
        """ Computes fracture leakage and bottom hole pressure at non-dimensional time t_d. 
        
        :param float t_d: Non-dimensional time. 
        :return (numpy.array, float): Returns a tuple, with the first entry \
        being the leakage over the length of the fracture. The \
        number of entries corresponds to the number_of_segments used for the
        semi-analytical solution. The second entry in the tuple is the bottom 
        hole pressure at that time. 
        """
        x_D = [float(i)/self.number_of_segments for i in range(self.number_of_segments)]

        delta_x = 1./self.number_of_segments

        def laplace_solution(s):
            LHS = np.zeros((self.number_of_segments+1, self.number_of_segments+1))
            RHS = np.zeros(self.number_of_segments+1)
            for j in range(self.number_of_segments):
                x_DJ = x_D[j]+.5*delta_x
                RHS[j] = np.pi*x_DJ/(self.c_fd*s)
                for i in range(self.number_of_segments):
                    x_DI = x_D[i]
                    if x_DJ > x_DI and x_DJ < x_DI+delta_x:
                        (quad1, error1) = integrate.quad(lambda x: special.kn(0, abs(x_DJ-x)*np.sqrt(s)), 
                                                         x_DI,
                                                         x_DI+delta_x,
                                                         points = [x_DJ], epsrel=1.e-16)
                    else:
                        (quad1, error1) = integrate.quad(lambda x: special.kn(0, abs(x_DJ-x)*np.sqrt(s)), 
                                                         x_DI,
                                                         x_DI+delta_x, epsrel=1.e-16)

                    (quad2, error2) = integrate.quad(lambda x: special.kn(0, abs(x_DJ+x)*np.sqrt(s)), 
                                                     x_DI,
                                                     x_DI+delta_x, epsrel=1.e-16)
                    #assert (error1 < 1.e-8)
                    #assert (error2 < 1.e-8)
                    LHS[j, i] += -.5*(quad1+quad2)
                    if i < j:
                        LHS[j, i] += (np.pi/self.c_fd)*((delta_x**2)/2.+
                                                   (x_DJ-(i+1)*delta_x)*delta_x)

                LHS[j, j] += (np.pi/self.c_fd)*(delta_x**2)/8.

            for i in range(self.number_of_segments):
                LHS[i, number_of_segments] += 1.

            for j in range(self.number_of_segments):
                LHS[number_of_segments, j] += delta_x

            RHS[self.number_of_segments] = 1./s
            solution = np.linalg.solve(LHS, RHS)
            min_singular_value = np.linalg.svd(LHS)[1][-1]
            assert(min_singular_value > 1.e-8)

            return solution

        full_solution = gaver_stehfest(t_D, laplace_solution)

        return (full_solution[:-1], full_solution[-1])

        
def cinco_ley_laplace(C_fD, t_D, number_of_segments = 20):
    """ Computes the semi-anlytical solution from 
    Cinco Ley and Meng:
    "Pressure Transient Analysis of Wells With Finite Conductivity 
    Vertical Fractures in Double Porosity Reservoirs"
    Also see:
    "Non-Darcy Flow in Wells With Finite-Conductivity 
    Vertical Fractures. 

    C_fD: Fracture conductivity difined as 
    (fracture_k*fracture_width)/(reservoir_k*fracture_length)

    t_D = \frac{\Beta res_k t}{\theta c_t \mu x_f^2}
    p_D = \frac{res_k h [p_i - p(t)]}{\alpha q \Beta \mu}
    q_D = \frac{2 q_f}{q_w}
    """
    x_D = [float(i)/number_of_segments for i in range(number_of_segments)]
    
    delta_x = 1./number_of_segments

    def laplace_solution(s):
        LHS = np.zeros((number_of_segments+1, number_of_segments+1))
        RHS = np.zeros(number_of_segments+1)
        for j in range(number_of_segments):
            x_DJ = x_D[j]+.5*delta_x
            RHS[j] = np.pi*x_DJ/(C_fD*s)
            for i in range(number_of_segments):
                x_DI = x_D[i]
                if x_DJ > x_DI and x_DJ < x_DI+delta_x:
                    (quad1, error1) = integrate.quad(lambda x: special.kn(0, abs(x_DJ-x)*np.sqrt(s)), 
                                                     x_DI,
                                                     x_DI+delta_x,
                                                     points = [x_DJ], epsrel=1.e-16)
                else:
                    (quad1, error1) = integrate.quad(lambda x: special.kn(0, abs(x_DJ-x)*np.sqrt(s)), 
                                                     x_DI,
                                                     x_DI+delta_x, epsrel=1.e-16)

                (quad2, error2) = integrate.quad(lambda x: special.kn(0, abs(x_DJ+x)*np.sqrt(s)), 
                                                 x_DI,
                                                 x_DI+delta_x, epsrel=1.e-16)
                #assert (error1 < 1.e-8)
                #assert (error2 < 1.e-8)
                LHS[j, i] += -.5*(quad1+quad2)
                if i < j:
                    LHS[j, i] += (np.pi/C_fD)*((delta_x**2)/2.+
                                               (x_DJ-(i+1)*delta_x)*delta_x)
                
            LHS[j, j] += (np.pi/C_fD)*(delta_x**2)/8.

        for i in range(number_of_segments):
            LHS[i, number_of_segments] += 1.

        for j in range(number_of_segments):
            LHS[number_of_segments, j] += delta_x
            
        RHS[number_of_segments] = 1./s
        solution = np.linalg.solve(LHS, RHS)
        min_singular_value = np.linalg.svd(LHS)[1][-1]
        assert(min_singular_value > 1.e-8)

        return solution
    
    full_solution = gaver_stehfest(t_D, laplace_solution)

    return (full_solution[:-1], full_solution[-1])
