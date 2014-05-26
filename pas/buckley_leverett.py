
import numpy as np
from scipy import sparse, diag
from scipy import interpolate

import matplotlib.pyplot as plt


class BuckleyLeverett():
    """ Solves the two-phase flow problem using the Buckley-Leverett 
    solution. The code provides three methods for providing the 
    specifying the fractional flow:
    1) df_ds (highest priority): Directly provide the derivative of the fractional flow \
    function with respect to saturation. 
    2) fractional_flow: Provide the fractinal flow function. 
    3) k_ro, k_rw, oil_viscosity, water_viscosity: Provide the variables and functions 
    that make up the fractional flow function. 

    The code first looks for a df_ds function, finding none is looks for fractional_flow, 
    and finally it will construct fractional_flow from the relative permeability and 
    fluid viscosities. 
    
    :ivar float cross_section: Domain cross sectional area (m^2). 
    :ivar float length: Domain length (m). 
    :ivar float injection_rate: Well water injection rate (m^3/s). 
    :ivar float residual_water: Residual water saturation. 
    :ivar float residual_oil: Residual oil saturation. 
    :ivar function fractional_flow: Fractional flow as a function \
    of water satuation (optional). 
    :ivar function df_ds_function: Derivative of fractional flow with \
    respect to saturation (optional). If not specified, the program \
    will compute a numerical derivative of fractional_flow. 
    """

    def __init__(self):
        
        ## Domain cross sectional area (m^2). 
        self.cross_section = 1.

        ## Domain length (m)
        self.length = 1.

        ## Rock porosity. 
        self.porosity = 1.
        
        ## Injection rate (m^3/s). 
        self.injection_rate = 1.

        ## Minimum water residual. 
        self.residual_water = 0.
        
        ## Minimum oil residual. 
        self.residual_oil = 0.
        
        ## Oil relative permeability
        self.k_ro = None

        ## Water relative permeability
        self.k_rw = None 
        
        ## Oil viscosity
        self.oil_viscosity = 1.e-3

        ## Water viscosity
        self.water_viscosity = 1.e-3

        ## Fractional flow as a function of water satuation ('function' type). 
        self.fractional_flow = None

        ## Derivative of fractional flow with respect to saturation ('function' type). 
        self.df_ds_function = None

    def plot_water_saturation_at_t(self, time):
        
        x = np.linspace(0., self.length, 50)
        y = [self.water_saturation_solution(time)(x_i) for x_i in x]
        
        plt.plot(x, y, 'k')
        plt.show()
        

    def water_saturation_solution(self, time):
        """ Returns a function based on the  Buckley-Leverett solution 
        problem for two-phase problems. The functions gives the water saturation 
        at the time (s) specified in the input. 
        :param float time: Time (s) after injection.
        :return: A function representing the water saturation at x (distance from injector in meters). 
        """
        if self.df_ds_function is None:
            if self.fractional_flow is None:
                def fractional_flow(sw):
                    return 1./(1.+self.k_ro(sw)/self.k_rw(sw)*self.oil_viscosity/self.water_viscosity)
                self.fractional_flow = fractional_flow
            def fractional_flow_prime(water_saturation):
                fw_prime = self.fractional_flow(water_saturation+.0001)
                fw_prime -=self.fractional_flow(water_saturation)
                fw_prime /= .0001
                return fw_prime
        else:
            fractional_flow_prime = self.df_ds_function

        (_, sw_max) = min(map(lambda sw:(abs(fractional_flow_prime(sw)-self.fractional_flow(sw)/sw), sw), 
                              np.arange(.01, 1.-.2, .0001)))

        f_prime_max = fractional_flow_prime(sw_max)

        v = time*self.injection_rate/(self.porosity*self.cross_section)
        x_front = v*fractional_flow_prime(sw_max)

        # Compute sw as a function of x. This is done by 
        # interpolating over sw, and finding f_prime, 
        # then scaling f_prime properly to cover x. 
        inter_y = [1.-self.residual_oil]+list(np.arange(1.-self.residual_oil, sw_max, -.001))+[sw_max]
        inter_x = [0.]
        inter_x += map(lambda sw: fractional_flow_prime(sw), 
                      np.arange(1.-self.residual_oil, sw_max, -.001))
        inter_x += [fractional_flow_prime(sw_max)]
        inter_x = map(lambda x:x/f_prime_max*x_front, inter_x)

        saturation_before_front = interpolate.interp1d(inter_x, inter_y)

        def bl_solution(point):
            if point<x_front:
                return saturation_before_front(point)
            else:
                return 0.

        return bl_solution

