
import pas.buckley_leverett as bl
import numpy as np

case1 = bl.BuckleyLeverett()

krw=lambda se:se**3
kro=lambda se:(1.-se)**3

case1.residual_n = .0
case1.residual_w = .0

krw = bl.rescale_relative_perm(krw, case1.residual_w, case1.residual_n)
kro = bl.rescale_relative_perm(kro, case1.residual_w, case1.residual_n)

case1.k_rw = krw
case1.k_rn = kro

case1.viscosity_w = 1.e-4
case1.viscosity_n = 1.e-4

case1.length = 5.

case1.plot_fractional_flow(show_sw_at_front = True)

(sol, front) = case1.saturation_solution(10.)

for x in np.arange(0., 5., .1):
    print x, sol(x)

case1.plot_saturation_at_t(1.)
case1.plot_saturation_at_t(2.)
#case1.plot_water_saturation_at_t(200.)
