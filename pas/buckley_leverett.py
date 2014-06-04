
import numpy as np
from scipy import sparse, diag
from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def rescale_relative_perm(k_r, residual_w, residual_n):
    """ Takes a relative permeability function from [0,1] and 
    rescales it to [residual_w, 1-residual_n]. 
    """
    def new_k_r(sw):
        se = sw-residual_w 
        se /= 1.-residual_w-residual_n
        return k_r(se)

    return new_k_r

class BuckleyLeverett():
    """ Solves the two-phase flow problem using the Buckley-Leverett 
    solution. 

    The code first looks for a df_ds function, finding none is looks 
    for fractional_flow, and finally it will construct fractional_flow 
    from the relative permeability and fluid viscosities. 
    
    :ivar float cross_section: Domain cross sectional area (m^2). 
    :ivar float length: Domain length (m). 
    :ivar float injection_rate: Injection rate (m^3/s). 
    :ivar float residual_w: Residual saturation for wetting phase. 
    :ivar float residual_n: Residual saturation for non-wetting phase. 
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

        ## Minimum wetting phase residual. 
        self.residual_w = 0.
        
        ## Minimum non-wetting phase residual. 
        self.residual_n = 0.
        
        ## Non-wetting phase relative permeability
        self.k_rn = None

        ## Wetting phase relative permeability
        self.k_rw = None 
        
        ## Non-wetting viscosity
        self.viscosity_n = 1.e-3

        ## Wetting viscosity
        self.viscosity_w = 1.e-3


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
                                            sw in np.arange(self.residual_w+1.e-12, 
                                                            1.-self.residual_n, 
                                                            dx)])
        df_ds_data = [(df_ds(sw), sw) for sw in np.arange(self.residual_w+1.e-12, 
                                                          1.-self.residual_n, 
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

    def plot_saturation_at_t(self, time):
        """
        """
        (solution, front)= self.saturation_solution(time)

        print "front", front
        
        if front < self.length:
            x = [x/100.*front for x in range(100)]+[front]
            y = [solution(x_i) for x_i in x]
            x += [front]
            y += [0.]
            x += [self.length]
            y += [0.]

        else:
            x = [x/100.*self.length for x in range(100)]+[self.length]
            y = [solution(x_i) for x_i in x]
            
        plt.plot(x, y, 'k')
        plt.show()

    def plot_fractional_flow(self, show_sw_at_front = False):
        """
        """
        fractional_flow = self.construct_fractional_flow()
        x = np.linspace(self.residual_w+1.e-3, 1.-self.residual_n, 100)
        y = [fractional_flow(x_i) for x_i in x]
        plt.plot(x, y, 'k')

        if show_sw_at_front:
            fractional_flow_prime = self.construct_fractional_flow_prime()
            sw_at_front = self.sw_at_front(fractional_flow, fractional_flow_prime)
            plt.plot([sw_at_front], [fractional_flow(sw_at_front)], "ro")
        
        plt.show()
        

    def plot_relative_permeability(self):
        """
        """
        sw_points = np.linspace(self.residual_w+1.e-3, 1.-self.residual_n, 100)
        k_rw_points = [self.k_rw(sw_i) for sw_i in sw_points]
        k_rn_points = [self.k_rn(sw_i) for sw_i in sw_points]
        

        plt.plot(sw_points, k_rw_points)
        plt.plot(sw_points, k_rn_points)

        plt.plot([1.-self.residual_n, 1.-self.residual_n], [0., self.k_rw(1.-self.residual_n)], 'k--')
        plt.plot([self.residual_w, self.residual_w], [0., self.k_rn(self.residual_w)],'k--' )

        plt.xlim([0., 1.])

        plt.show()
        
    def plot_fractional_flow_prime(self):
        """
        """
        fractional_flow_prime = self.construct_fractional_flow_prime()
        x = np.linspace(self.residual_w+1.e-3, 1.-self.residual_n, 100)
        y = [fractional_flow_prime(x_i) for x_i in x]
        plt.plot(x, y, 'k')
        plt.show()

    def sw_at_front(self, fractional_flow, fractional_flow_prime):
        """ Find the saturation at the front. 
        """
        
        fractional_flow_prime_prime = self.construct_fractional_flow_prime_prime()
        
        ## Find the sw range to search in. This is done using the region 
        ## the second derivative is negative. 
        sw_start = 1.-self.residual_n 
        sw_end = self.residual_w

        found_start = False
        is_linear = True
        
        for sw in np.arange(self.residual_w+1.e-10, 1.-self.residual_n, .001):
            if abs(fractional_flow_prime_prime(sw)) > 1.e-3:
                is_linear = False
            if fractional_flow_prime_prime(sw) < -1.e-2 and sw < sw_start: 
                sw_start = sw
            if fractional_flow_prime_prime(sw) < -1.e-2 and sw > sw_end: 
                sw_end = sw

        if is_linear:
            sw_at_front = 1.-self.residual_n

        else:            
            (_, sw_at_front) = min(map(lambda sw:(abs(fractional_flow_prime(sw)-fractional_flow(sw)/sw), sw), 
                                       np.arange(sw_start, 
                                                 sw_end, 
                                                 .0001)))
        
        return sw_at_front

    def construct_fractional_flow(self):
        """ Builds fractional flow function from 
        parameters of problem. 
        """
        def fractional_flow(sw):
            return 1./(1.+self.k_rn(sw)/self.k_rw(sw)*\
                           self.viscosity_n/self.viscosity_w)
        return fractional_flow

    def construct_fractional_flow_prime(self):
        """ Builds fractional flow  derivative function from 
        parameters of problem. 
        """
        fractional_flow = self.construct_fractional_flow()
        
        def fractional_flow_prime(sw):
            f_prime = fractional_flow(sw+.0001)
            f_prime -=fractional_flow(sw)
            f_prime /= .0001
            return f_prime
    
        return fractional_flow_prime

    def construct_fractional_flow_prime_prime(self):
        """ Builds fractional flow  derivative function from 
        parameters of problem. 
        """
        fractional_flow = self.construct_fractional_flow()
        
        def fractional_flow_prime_prime(sw):
            h = .0001
            f_prime = fractional_flow(sw+h)
            f_prime -=2.*fractional_flow(sw)
            f_prime += fractional_flow(sw-h)
            f_prime /= h**2
            return f_prime
        
        return fractional_flow_prime_prime

    def saturation_solution(self, time):
        """ Returns the Buckley-Leverett solution 
        for two-phase problems. The functions gives the wetting phase 
        saturation at the time (s) specified in the input. 
        :param float time: Time (s) after injection.
        :return: A function representing the wetting phase saturation at x \
        (distance from injector in meters). 
        """
        fractional_flow = self.construct_fractional_flow()
        fractional_flow_prime = self.construct_fractional_flow_prime()

        sw_at_front = self.sw_at_front(fractional_flow, fractional_flow_prime)

        fractional_flow_prime_at_front = fractional_flow_prime(sw_at_front)

        v = time*self.injection_rate/(self.porosity*self.cross_section)
        
        x_front = v*fractional_flow_prime(sw_at_front)

        # Compute sw as a function of x. This is done by 
        # interpolating over sw, and finding f_prime, 
        # then scaling f_prime properly to cover x. 
        inter_y = [1.-self.residual_n]+list(np.arange(1.-self.residual_n, sw_at_front, -.001))+[sw_at_front]
        inter_x = [0.]
        inter_x += map(lambda sw: fractional_flow_prime(sw), 
                      np.arange(1.-self.residual_n, sw_at_front, -.001))
        inter_x += [fractional_flow_prime(sw_at_front)]
        inter_x = map(lambda x:x/fractional_flow_prime_at_front*x_front, inter_x)

        saturation_before_front = interpolate.interp1d(inter_x, inter_y)

        def bl_solution(point):
            if point<=x_front:
                return saturation_before_front(point)
            else:
                return 0.

        return (bl_solution, x_front)

