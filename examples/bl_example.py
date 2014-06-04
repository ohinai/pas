
import pas.buckley_leverett as bl
import numpy as np

case1 = bl.BuckleyLeverett()

krw=lambda se:.4*se**2
kro=lambda se:(1.-se)**2

case1.residual_n = .2
case1.residual_w = .2

krw = bl.rescale_relative_perm(krw, case1.residual_w, case1.residual_n)
kro = bl.rescale_relative_perm(kro, case1.residual_w, case1.residual_n)

case1.k_rw = krw
case1.k_rn = kro

case1.viscosity_w = 1.e-4
case1.viscosity_n = 1.e-4

case1.length = 5.

case1.plot_relative_permeability()

#case1.plot_fractional_flow(show_sw_at_front = True)
#case1.plot_fractional_flow_prime()

(sol, front) = case1.saturation_solution(10.)

case1.plot_saturation_at_t(1.)
#case1.plot_water_saturation_at_t(200.)
