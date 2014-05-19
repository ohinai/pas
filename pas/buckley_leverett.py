
import numpy as np
from scipy import sparse, diag
from scipy import interpolate

def buckley_leverett_solution(self, current_time):
    """ Returns a function with solution to the Buckley-Leverett 
    problem for special two-phase problems. 
    """
    def fractional_flow(water_saturation):
        return self.water_mobility(water_saturation, self.ref_pressure_oil)/\
            (self.oil_mobility(water_saturation, self.ref_pressure_oil)+
             self.water_mobility(water_saturation, self.ref_pressure_oil))

    def fractional_flow_prime(water_saturation):
        fw_prime = fractional_flow(water_saturation+.000001)
        fw_prime -=fractional_flow(water_saturation)
        fw_prime /= .000001
        return fw_prime

    (_, sw_max) = min(map(lambda sw:(abs(fractional_flow_prime(sw)-fractional_flow(sw)/sw), sw), 
                          np.arange(.01, 1.-.2, .0001)))

    # The code assumes the first well to be the only well 
    # placed on a side of the domain. 
    Q_t = self.rate_wells_rate_water[0]/self.ref_density_water
    f_prime_max = fractional_flow_prime(sw_max)

    # Assume the height is in the z-direction. 
    A = self.mesh.get_dim_z()
    v = current_time*Q_t/(self.porosities[0]*A)
    x_front = v*fractional_flow_prime(sw_max)

    # Figure out sw as a function of x. This is done by 
    # interpolating over sw, and finding f_prime, 
    # then scaling f_prime properly to cover x. 
    inter_y = [1.-self.residual_saturation_oil]+list(np.arange(1.-self.residual_saturation_oil, sw_max, -.001)) + [sw_max]
    inter_x = [0.] 
    inter_x += map(lambda sw: fractional_flow_prime(sw), 
                  np.arange(1.-self.residual_saturation_oil, sw_max, -.001))
    inter_x += [fractional_flow_prime(sw_max)]
    inter_x = map(lambda x:x/f_prime_max*x_front, inter_x)

    saturation_before_front = interpolate.interp1d(inter_x, inter_y)

    def bl_solution(point):
        if point[0]<x_front:
            return saturation_before_front(point[0])
        else:
            return 0.

    return bl_solution
