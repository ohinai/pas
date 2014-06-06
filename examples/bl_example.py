
import pas.buckley_leverett as bl
import numpy as np

case1 = bl.BuckleyLeverett()

case1.parameter["residual_n"] = .2
case1.parameter["residual_w"] = .1
case1.parameter["initial_sw"] = .1

def krw(sw):
    se = sw-case1.parameter["residual_w"]
    se /= 1.-case1.parameter["residual_n"]-case1.parameter["residual_w"]
    return .4*se**4

def krn(sw):
    se = sw-case1.parameter["residual_w"]
    se /= 1.-case1.parameter["residual_n"]-case1.parameter["residual_w"]
    return (1.-se)**4

case1.k_rw = krw
case1.k_rn = krn

case1.parameter["viscosity_w"] = 1.e-4
case1.parameter["viscosity_n"] = 1.e-4

case1.parameter["length"] = 5.

case1.plot_saturation_at_t([1., 2., 10.])

case1.plot_relative_permeability()
case1.plot_fractional_flow(show_sw_at_front = True)
