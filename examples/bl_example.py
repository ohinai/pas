
import pas.buckley_leverett as bl
import numpy as np

case1 = bl.BuckleyLeverett()

krw=lambda se:se**1
kro=lambda se:(1.-se)**1

case1.residual_oil = .0
case1.residual_water = .0

krw = bl.rescale_relative_perm(krw, case1.residual_water, case1.residual_oil)
kro = bl.rescale_relative_perm(kro, case1.residual_water, case1.residual_oil)

water_viscosity = 1.e-4
oil_viscosity = 1.e-4

def fractional_flow(water_saturation):
    ff_value = krw(water_saturation)/water_viscosity
    ff_value /= ff_value+kro(water_saturation)/oil_viscosity
    return ff_value

case1.fractional_flow = fractional_flow
case1.length = 5.

case1.plot_fractional_flow()

sol = case1.water_saturation_solution(10.)

for x in np.arange(0., 5., .1):
    print x, sol(x)

case1.plot_water_saturation_at_t(1.)
#case1.plot_water_saturation_at_t(200.)
