
import numpy as np
from scipy import sparse, diag
from scipy import interpolate

class BuckleyLeverett():
    """ Solves the two-phase flow problem using the Buckley-Leverett 
    solution. 
    
    :ivar float height: Domain height (m). 
    :ivar float injection_rate: Well water injection rate (kg/s). 
    :ivar float residual_water: Residual water saturation. 
    :ivar float residual_oil: Residual oil saturation. 
    :ivar function fractional_flow: Fractional flow as a function \
    of water satuation. 
    :ivar function df_ds_function: Derivative of fractional flow with \
    respect to saturation (optional). If not specified, the program \
    will compute a numerical derivative of fractional_flow. 
    """

    def __init__(self):
        
        ## Domain height (m). 
        self.height = 1.

        ## Rock porosity. 
        self.porosity = 1.
        
        ## Injection rate (kg/s). 
        self.injection_rate = 1.

        ## Minimum water residual. 
        self.residual_water = 0.
        
        ## Minimum oil residual. 
        self.residual_oil = 0.
        
        ## Fractional flow as a function of water satuation ('function' type). 
        self.fractional_flow = None

        ## Derivative of fractional flow with respect to saturation ('function' type). 
        self.df_ds_function = None

    def water_saturation_solution(self, time):
        """ Returns a function based on the  Buckley-Leverett solution 
        problem for two-phase problems. The functions gives the water saturation 
        at the time (s) specified in the input. 
        :param float time: Time (s) after injection.
        :return: A function representing the water saturation at x (distance from injector in meters). 
        """
        
        if self.df_ds_function is None:
            def fractional_flow_prime(water_saturation):
                fw_prime = self.fractional_flow(water_saturation+.000001)
                fw_prime -=self.fractional_flow(water_saturation)
                fw_prime /= .000001
                return fw_prime
        else:
            fractional_flow_prime = self.df_ds_function

        (_, sw_max) = min(map(lambda sw:(abs(fractional_flow_prime(sw)-self.fractional_flow(sw)/sw), sw), 
                              np.arange(.01, 1.-.2, .0001)))

        f_prime_max = fractional_flow_prime(sw_max)

        v = time*self.injection_rate/(self.porosity*self.height)
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
