
import numpy as np
from scipy import sparse, diag
from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def rescale_relative_perm(k_r, residual_water, residual_oil):
    """ Takes a relative permeability function from [0,1] and 
    rescales it to [residual_water, 1-residual_oil]. 
    """
    def new_k_r(sw):
        se = sw-residual_water 
        se /= 1.-residual_water-residual_oil
        return k_r(se)

    return new_k_r

class BuckleyLeverett():
    """ Solves the two-phase flow problem using the Buckley-Leverett 
    solution. The code provides three methods for providing the 
    specifying the fractional flow:
    1) df_ds (highest priority): Directly provide the derivative of \
    the fractional flow \
    function with respect to saturation. 
    2) fractional_flow: Provide the fractinal flow function. 
    3) k_ro, k_rw, oil_viscosity, water_viscosity: Provide the \
    variables and functions 
    that make up the fractional flow function. 

    The code first looks for a df_ds function, finding none is looks 
    for fractional_flow, and finally it will construct fractional_flow 
    from the relative permeability and fluid viscosities. 
    
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

    

    def find_df_ds_max_by_intergration(self, df_ds):
        """ Find df_ds max using the original method 
        proposed by Buckley-Leverett. 

        :param function df_ds: Function representing the derivative 
        of the fractional flow with respect to saturation. 
        :return: 
        """
        ## Find the half-way point on the curve, the point 
        ## at which is starts decreasing. This point 
        ## is located by finding the max value of df_ds. 

        dx = .0001
        (fw_mid, sw_mid) = max([(df_ds(sw), sw) for \
                                            sw in np.arange(self.residual_water+1.e-12, 
                                                            1.-self.residual_oil, 
                                                            dx)])
        df_ds_data = [(df_ds(sw), sw) for sw in np.arange(self.residual_water+1.e-12, 
                                                          1.-self.residual_oil, 
                                                          dx)]
                                          
        def find_difference_in_volume(f_prime_max):
            integral_1 = 0.
            integral_2 = 0.
            sw_max = 0.
            for (f_p, sw) in df_ds_data:
                if f_p < f_prime_max and sw < sw_mid:
                    integral_1 += f_p
                    sw_max = sw
                elif f_p > f_prime_max:
                    integral_2 += f_p
                    sw_max_2 = sw

            integral_1 = f_prime_max*sw_max - integral_1*dx
            integral_2 *= dx
            integral_2 -= f_prime_max*(sw_max_2-sw_max)
            
            return (abs(integral_1-integral_2), sw_max, sw_max_2)

        current_min = 1000.
        current_sw = 0.
        plt.plot([x[1] for x in df_ds_data], [x[0] for x in df_ds_data], 'k')
        plt.show()
        for (f_p, sw) in df_ds_data:
            if sw >sw_mid:
                (diff, current_sw1, current_sw2) = find_difference_in_volume(f_p)
                if diff < current_min:
                    current_sw = sw
                    current_min = diff
                    sw1 = current_sw1
                    sw2 = current_sw2
                    
        print current_sw

        fig, ax = plt.subplots()

        plt.plot([x[1] for x in df_ds_data], [x[0] for x in df_ds_data], 'k')
        plt.plot([sw_mid, current_sw], [fw_mid, df_ds(current_sw)], 'ro')
        plt.plot([sw1], [df_ds(sw1)], 'ro')
        
        verts1 = [ (x[1], x[0]) for x in filter(lambda x:x[1]<sw1, df_ds_data)]
        verts1 += [(sw1, df_ds(sw1)), (0., df_ds(current_sw)), (0.,0.)]
        fill1 = Polygon(verts1, facecolor='0.9', edgecolor='0.5')
        ax.add_patch(fill1)

        verts2 = [ (x[1], x[0]) for x in filter(lambda x:x[1]>sw1 and x[1]<sw2, df_ds_data)]
        fill2 = Polygon(verts2, facecolor='0.9', edgecolor='0.5')
        ax.add_patch(fill2)

        plt.grid(True)
        plt.show()

    def plot_water_saturation_at_t(self, time):
        
        x = np.linspace(0., self.length, 100)
        solution = self.water_saturation_solution(time)
        y = [solution(x_i) for x_i in x]
        
        plt.plot(x, y, 'k')
        plt.show()

    def plot_fractional_flow(self):

        x = np.linspace(self.residual_water, 1.-self.residual_oil, 100)
        y = [self.fractional_flow(x_i) for x_i in x]
        
        plt.plot(x, y, 'k')
        plt.show()
        

    def water_saturation_solution(self, time):
        """ Returns a function based on the  Buckley-Leverett solution 
        for two-phase problems. The functions gives the water saturation 
        at the time (s) specified in the input. 
        :param float time: Time (s) after injection.
        :return: A function representing the water saturation at x \
        (distance from injector in meters). 
        """
        if self.df_ds_function is None:
            if self.fractional_flow is None:
                def fractional_flow(sw):
                    return 1./(1.+self.k_ro(sw)/self.k_rw(sw)*\
                                   self.oil_viscosity/self.water_viscosity)
                self.fractional_flow = fractional_flow
            def fractional_flow_prime(water_saturation):
                fw_prime = self.fractional_flow(water_saturation+.0001)
                fw_prime -=self.fractional_flow(water_saturation)
                fw_prime /= .0001
                return fw_prime
        else:
            fractional_flow_prime = self.df_ds_function

        (_, sw_max) = min(map(lambda sw:(abs(fractional_flow_prime(sw)-self.fractional_flow(sw)/sw), sw), 
                              np.arange(self.residual_water+.01, 
                                        1.-self.residual_oil, 
                                        .0001)))
        print "sw_max", sw_max

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

