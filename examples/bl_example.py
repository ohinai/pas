
import pas.buckley_leverett as bl
import numpy as np

case1 = bl.BuckleyLeverett()

krw=lambda se:se**2 
kro=lambda se:(1.-se)**2

water_viscosity = 1.e-4
oil_viscosity = 1.e-4

def fractional_flow(water_saturation):
    ff_value = krw(water_saturation)/water_viscosity
    ff_value /= ff_value+kro(water_saturation)/oil_viscosity
    return ff_value

case1.fractional_flow = fractional_flow
sol = case1.water_saturation_solution(2.)

for x in np.arange(0, 5., .1):
    print x, sol(x)
