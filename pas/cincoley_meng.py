

import numpy as np
from scipy import sparse, diag, integrate, special 
import math as math
import scipy.sparse.linalg.dsolve as dsolve
import os
from multiprocessing import Pool
import itertools

def cinco_ley_solution_time():
    """ Computes the semi-anlytical solution from 
    Cinco Ley et al. 
    """
    # The construction comes from the original 
    # paper by Heber Cinco L. et al 
    # "Transient Pressure Behavior for the Well 
    # With a Finite-Conductivity Vertical Fracture."
    # The notation resembles the variables used 
    # in equation (11). 
    # special.expi => E_i
    # math.erf => erf

    # The paper has the numbering start from 1.
    # Fracture segment indices
    N = 10

    # Number of time segments (including time 0)
    K = 2

    delta_t = .1
    delta_x = .1 
    
    w = .001

    intial_p = 100.

    frac_phi = 1.
    res_phi = .2

    compressibility = 1.e-9
    
    res_k = 1.e-12
    frac_k = 1.e-5
    
    x_f = 10.
    
    C_fDf = w*frac_phi*compressibility
    C_fDf /= np.pi*x_f*res_phi*compressibility
    
    eta_fD = frac_k*frac_phi*compressibility
    eta_fD /= res_k*res_phi*compressibility

    t_D = [(i)*delta_t for i in range(K)]
    x_D = [(i+.5)/N for i in range(N)]
    
    RHS = np.zeros(N)
    LHS = np.zeros((N, N))

    def alpha(i, j):
        return (j-i+2.+.5)/(2.*N)

    def beta(i, j):
        return (j-i+2.-.5)/(2.*N)
    
    def gamma(i, j):
        # Not sure if it's i or 1. 
        # There might be an error 
        # in the paper. 
#        return (j+i+2.-.5)/(2.*N)
        return (j+1.+1.-.5)/(2.*N)

    def delta(i, j):
        return (j+i+2.-1.5)/(2.*N)

    def X(l, i, j):
        x_value = math.erf(alpha(i, j)/(np.sqrt(t_D[K-1]-t_D[l])))
        x_value -= math.erf(beta(i, j)/(np.sqrt(t_D[K-1]-t_D[l])))
        x_value += math.erf(gamma(i,j)/(np.sqrt(t_D[K-1]-t_D[l])))
        x_value -= math.erf(delta(i,j)/(np.sqrt(t_D[K-1]-t_D[l])))
        x_value *= 2.*np.sqrt(t_D[K-1]-t_D[l])
        return x_value

    def Y(l, i, j):
        y_value = alpha(i, j)*special.expi(-alpha(i, j)**2/(t_D[K-1]-t_D[l]))
        y_value -= beta(i, j)*special.expi(-beta(i, j)**2/(t_D[K-1]-t_D[l]))
        y_value += gamma(i, j)*special.expi(-gamma(i, j)**2/(t_D[K-1]-t_D[l]))
        y_value -= delta(i, j)*special.expi(-delta(i, j)**2/(t_D[K-1]-t_D[l]))
        y_value *= -2./np.sqrt(np.pi)
        return y_value
    
    l = 1
    print "C_fDf", C_fDf
    print "eta_fD", eta_fD
    
    for j in range(N):
        RHS[j] = t_D[K-1]
        RHS[j] += 2./(np.pi**2*eta_fD)*analytical.cincoley_sum_1(.001, x_D[j], eta_fD, 100000)
        # sum([1./n**2*(1.-np.exp(-np.pi**2*n*n*eta_fD*t_D[K-1]))*\
        #          np.cos(n*np.pi*x_D[j])\
        #          for n in xrange(1, 1000000)])
        RHS[j] *= 1./C_fDf

        for i in range(N):
            # print sum([1./n**3*\
            #                (np.exp(-np.pi**2*n*n*eta_fD*(t_D[K-1]-t_D[l]))-\
            #                     np.exp(-np.pi**2*n*n*eta_fD*(t_D[K-1]-t_D[l-1])))*\
            #                np.cos(n*np.pi*x_D[j])*np.cos(n*np.pi*x_D[i])*np.sin(n*np.pi/(2.*N))\
            #                for n in range(1, 10000)]), 
            LHS[j, i] = -4./(np.pi**3*eta_fD)*\
                analytical.cincoley_sum_2(t_D[K-1]-t_D[l], t_D[K-1]-t_D[l-1],\
                                              x_D[j], x_D[i], eta_fD, N,  100000)

            LHS[j, i] -= delta_x*t_D[l]-t_D[l-1]
            LHS[j, i] *= 1./C_fDf
            LHS[j, i] -= np.sqrt(np.pi)/4.*(X(l-1, i, j)-X(l, i, j)+\
                Y(l-1, i, j)-Y(l, i, j))

    print LHS
    print RHS 
    print "singular values", np.linalg.svd(LHS)[1]
    print np.linalg.solve(LHS, RHS)


def stehfest_laplace_inverse(t, P):
    """ Performs a numerical inversion 
    of the Laplace transform based 
    on the Stehfest algorithm 
    found in Well Test Analysis 
    by Rajagopal Raghavan on p. 90. 
    """
    N = 10
    def v_i(i):
        summation = 0
        for k in range(int((i+1.)/2), min(i, N/2)):
            current_sum = pow(k, N/2)*math.factorial(2*k)
            current_sum /= math.factorial(N/2-k)
            current_sum /= math.factorial(k)
            current_sum /= math.factorial(k-1)
            current_sum /= math.factorial(i-k)
            current_sum /= math.factorial(2*k-i)
            
            summation += current_sum
        return pow(-1, N/2+i)

    total_sum = v_i(1)*P(1.*np.log(2)/t)
    for i in range(2, 50):
        total_sum += v_i(i)*P(i*np.log(2.)/t)
    return total_sum*np.log(2.)/t

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
