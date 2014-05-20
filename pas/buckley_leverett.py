
import numpy as np
from scipy import sparse, diag
from scipy import interpolate

def buckley_leverett(current_time, 
                     fractional_flow, 
                     residual_water, 
                     residual_oil,
                     A, 
                     injection_rate, 
                     porosity):
    """ Returns a function with solution to the Buckley-Leverett 
    problem for two-phase problems. 
    """
    def fractional_flow_prime(water_saturation):
        fw_prime = fractional_flow(water_saturation+.000001)
        fw_prime -=fractional_flow(water_saturation)
        fw_prime /= .000001
        return fw_prime

    (_, sw_max) = min(map(lambda sw:(abs(fractional_flow_prime(sw)-fractional_flow(sw)/sw), sw), 
                          np.arange(.01, 1.-.2, .0001)))

    Q_t = injection_rate
    f_prime_max = fractional_flow_prime(sw_max)

    v = current_time*Q_t/(porosity*A)
    x_front = v*fractional_flow_prime(sw_max)

    # Compute sw as a function of x. This is done by 
    # interpolating over sw, and finding f_prime, 
    # then scaling f_prime properly to cover x. 
    inter_y = [1.-residual_oil]+list(np.arange(1.-residual_oil, sw_max, -.001))+[sw_max]
    inter_x = [0.]
    inter_x += map(lambda sw: fractional_flow_prime(sw), 
                  np.arange(1.-residual_oil, sw_max, -.001))
    inter_x += [fractional_flow_prime(sw_max)]
    inter_x = map(lambda x:x/f_prime_max*x_front, inter_x)

    saturation_before_front = interpolate.interp1d(inter_x, inter_y)

    def bl_solution(point):
        if point<x_front:
            return saturation_before_front(point)
        else:
            return 0.

    return bl_solution
