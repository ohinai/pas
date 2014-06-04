
import pas.buckley_leverett as bl
import numpy as np

case1 = bl.BuckleyLeverett()

case1.residual_n = .2
case1.residual_w = .0
case1.initial_sw = .0

def krw(sw):
    se = sw-case1.residual_w
    se /= 1.-case1.residual_n-case1.residual_w
    return se**1

def krn(sw):
    se = sw-case1.residual_w
    se /= 1.-case1.residual_n-case1.residual_w
    return (1.-se)**1


case1.k_rw = krw
case1.k_rn = krn

case1.viscosity_w = 1.e-4
case1.viscosity_n = 1.e-4

case1.length = 5.

#case1.plot_relative_permeability()

#case1.plot_fractional_flow(show_sw_at_front = True)
#case1.plot_fractional_flow_prime()

(sol, front) = case1.saturation_solution(10.)

case1.plot_saturation_at_t(1.)
#case1.plot_water_saturation_at_t(200.)
