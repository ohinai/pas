""" Module for solving and visualizing the Buckley-Leverett solution
for two-phase flow problems.
"""
import numpy as np
from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

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

        self.param = {
            ## Domain cross sectional area (m^2).
            "cross_section":1.,
            ## Domain length (m).
            "length":1.,
            ## Rock porosity.
            "porosity":1.,
            ## Injection rate (m^3/s).
            "injection_rate":1.,
            ## Minimum wetting phase residual.
            "residual_w":0.,
            ## Minimum non-wetting phase residual.
            "residual_n":0.,
            ## Initial saturation.
            "initial_sw":0.,
            ## Non-wetting viscosity.
            "viscosity_n":1.e-3,
            ## Wetting viscosity.
            "viscosity_w":1.e-3}

    def k_rn(self, sw):
        """ Relative permeability for non-wetting phase.
        To be defined by user.
        """
        raise NotImplementedError("k_rn function not defined.")

    def k_rw(self, sw):
        """ Relative permeability for wetting phase.
        To be defined by user.
        """
        raise NotImplementedError("k_rw function not defined.")

    def sw_at_front_intergration(self):
        """ Alternative method for finding the saturation
        at the front using integration of the derivative
        of the fractional flow function.

        :return:
        """
        ## Find the half-way point on the curve, the point
        ## at which is starts decreasing. This point
        ## is located by finding the max value of df_ds.
        dx = .0001
        fw_mid = 0.
        sw_mid = 0.
        ff_p_min = 1000.
        for sw in np.arange(self.param["residual_w"]+1.e-12,
                            1.-self.param["residual_n"],
                            dx):
            if self.fractional_flow_p(sw) < ff_p_min:
                ff_p_min = fractional_flow_p(sw)
                sw_mid = sw
                fw_mid = self.fractional_flow_p(sw)

        df_ds_data = [(self.fractional_flow_p(sw), sw) \
                          for sw in np.arange(self.param["residual_w"]+1.e-12,
                                              1.-self.param["residual_n"],
                                              dx)]

        def areas_at_f_prime(f_prime_max):
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
        sw_at_front = 0.

        for (f_p, sw) in df_ds_data:
            if sw > sw_mid:
                (diff, current_sw1, current_sw2) = areas_at_f_prime(f_p)
                if diff < current_min:
                    sw_at_front = sw
                    current_min = diff
                    sw1 = current_sw1
                    sw2 = current_sw2

        fig, ax = plt.subplots()

        plt.plot([x[1] for x in df_ds_data], [x[0] for x in df_ds_data], 'k')
        plt.plot([sw_mid], [fw_mid], 'ro')
        plt.plot([sw_at_front], [self.fractional_flow_p(sw_at_front)], 'ro')
        plt.plot([sw1], [self.fractional_flow_p(sw1)], 'ro')

        verts1 = []
        for (f_p, sw) in df_ds_data:
            if sw < sw1:
                verts1.append((sw, f_p))

        verts1 += [(sw1, self.fractional_flow_p(sw1)),
                   (0., self.fractional_flow_p(sw_at_front)), (0.,0.)]
        fill1 = Polygon(verts1, facecolor='0.9', edgecolor='0.5')
        ax.add_patch(fill1)

        verts2 = []
        for (f_p, sw) in df_ds_data:
            if sw > sw1 and sw < sw2:
                verts2.append((sw, f_p))

        fill2 = Polygon(verts2, facecolor='0.9', edgecolor='0.5')
        ax.add_patch(fill2)

        plt.grid(True)
        plt.show()

        return sw_at_front

    def plot_saturation_at_t(self, time_list):
        """ Plots the Buckley-Leverett solution at the different times
        listed in time_list.

        :param list time_list: List of times to solve for. The different
        plots will be plotted in a single graph.
        :return: None
        """
        for time in time_list:
            (solution, front)= self.saturation_solution(time)
            if front < self.param["length"]:
                x_axis = [x/100.*front for x in range(100)]+[front]
                y_axis = [solution(x_i) for x_i in x_axis]
                x_axis += [front]
                y_axis += [self.param["initial_sw"]]
                x_axis += [self.param["length"]]
                y_axis += [self.param["initial_sw"]]

            else:
                x_axis = [x/100.*self.param["length"] for x in range(100)]
                x_axis += [self.param["length"]]
                y_axis = [solution(x_i) for x_i in x_axis]

            plt.plot(x_axis, y_axis, 'k')

        plt.show()

    def plot_fractional_flow(self, show_sw_at_front = False):
        """ Plots the fractional flow curve based on the
        problem parameters given.

        :param bool show_sw_at_front: Draws a circle at the
        saturation calculated at the front.
        :return: None
        """
        x_axis = np.linspace(self.param["residual_w"]+1.e-3,
                             1.-self.param["residual_n"], 100)
        y_axis = [self.fractional_flow(x_i) for x_i in x_axis]
        plt.plot(x_axis, y_axis, 'k')

        if show_sw_at_front:
            sw_at_front = self.sw_at_front()
            plt.plot([sw_at_front], [self.fractional_flow(sw_at_front)], "ro")

        plt.show()

    def plot_relative_permeability(self):
        """ Plots the relative permeability curves for the wetting
        and non-wetting phases.
        """
        sw_points = np.linspace(self.param["residual_w"]+1.e-3,
                                1.-self.param["residual_n"], 100)
        k_rw_points = [self.k_rw(sw_i) for sw_i in sw_points]
        k_rn_points = [self.k_rn(sw_i) for sw_i in sw_points]

        plt.plot(sw_points, k_rw_points)
        plt.plot(sw_points, k_rn_points)

        plt.plot([1.-self.param["residual_n"],
                  1.-self.param["residual_n"]], \
                     [0., self.k_rw(1.-self.param["residual_n"])], 'k--')
        plt.plot([self.param["residual_w"], self.param["residual_w"]], \
                     [0., self.k_rn(self.param["residual_w"])],'k--' )

        plt.xlim([0., 1.])

        plt.show()

    def plot_fractional_flow_p(self):
        """ Plots the derivative of fractional flow
        with respect to saturation (df_ds).
        """
        x_axis = np.linspace(self.param["residual_w"]+1.e-3,
                             1.-self.param["residual_n"], 100)
        y_axis = [self.fractional_flow_p(x_i) for x_i in x_axis]
        plt.plot(x_axis, y_axis, 'k')
        plt.show()

    def sw_at_front(self):
        """ Find the saturation at the front.
        """
        ## Find the sw range to search in. This is done using the region
        ## the second derivative is negative.
        sw_start = 1.-self.param["residual_n"]
        sw_end = self.param["residual_w"]

        is_linear = True

        for sw in np.arange(self.param["residual_w"]+1.e-10,
                            1.-self.param["residual_n"], .001):
            if abs(self.fractional_flow_p_p(sw)) > 1.e-3:
                is_linear = False
            if self.fractional_flow_p_p(sw) < -1.e-2 and sw < sw_start:
                sw_start = sw
            if self.fractional_flow_p_p(sw) < -1.e-2 and sw > sw_end:
                sw_end = sw

        if is_linear:
            sw_at_front = 1.-self.param["residual_n"]

        else:
            sw_at_front = 0.
            current_min = 1000.
            for sw in np.arange(sw_start, sw_end, .0001):
                current_diff = abs(self.fractional_flow_p(sw)-\
                                       self.fractional_flow(sw)/sw)
                if current_diff < current_min:
                    sw_at_front = sw
                    current_min = current_diff

        return sw_at_front

    def fractional_flow(self, sw):
        """ Returns fractional flow function from
        parameters of problem.

        :param float sw: Saturation of wetting phase.
        :return: float
        """
        return 1./(1.+self.k_rn(sw)/self.k_rw(sw)*\
                       self.param["viscosity_n"]/self.param["viscosity_w"])

    def fractional_flow_p(self, sw):
        """ Returns the fractional flow derivative function from
        parameters of problem.

        :param float sw: Saturation of wetting phase.
        :return: float
        """
        f_prime = self.fractional_flow(sw+.0001)
        f_prime -= self.fractional_flow(sw)
        f_prime /= .0001
        return f_prime

    def fractional_flow_p_p(self, sw):
        """ Returns the fractional flow  second derivative function from
        parameters of problem.

        :param float sw: Saturation of wetting phase.
        :return: float
        """
        h = .0001
        f_prime = self.fractional_flow(sw+h)
        f_prime -= 2.*self.fractional_flow(sw)
        f_prime += self.fractional_flow(sw-h)
        f_prime /= h**2
        return f_prime

    def saturation_solution(self, time):
        """ Returns the Buckley-Leverett solution
        for two-phase problems. The functions gives the wetting phase
        saturation at the time (s) specified in the input.
        :param float time: Time (s) after injection.
        :return: A function representing the wetting phase saturation at x \
        (distance from injector in meters).
        """
        sw_at_front = self.sw_at_front()

        fractional_flow_p_at_front = self.fractional_flow_p(sw_at_front)

        front_v = time
        front_v *= self.param["injection_rate"]
        front_v /= self.param["porosity"]*self.param["cross_section"]

        x_front = front_v*self.fractional_flow_p(sw_at_front)

        # Compute sw as a function of x. This is done by
        # interpolating over sw, and finding f_prime,
        # then scaling f_prime properly to cover x.
        inter_y = [1.-self.param["residual_n"]]
        inter_y += np.arange(1.-self.param["residual_n"], sw_at_front, -.001)
        inter_y += [sw_at_front]
        inter_x = [0.]
        inter_x += [self.fractional_flow_p(sw) for sw in \
                      np.arange(1.-self.param["residual_n"], 
                                sw_at_front, -.001)]
        inter_x += [self.fractional_flow_p(sw_at_front)]
        inter_x = [x/fractional_flow_p_at_front*x_front for x in  inter_x]

        saturation_before_front = interpolate.interp1d(inter_x, inter_y)

        def bl_solution(point):
            if point<=x_front:
                return saturation_before_front(point)
            else:
                return self.param["initial_sw"]

        return (bl_solution, x_front)

