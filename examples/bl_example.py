
import os
pas_absolute_path = os.path.abspath("../")
import sys 
sys.path.append(pas_absolute_path)

import pas.buckley_leverett as bl

krw=lambda se:se 
kro=lambda se:(1.-se)

water_viscosity = 1.e-4
oil_viscosity = 1.e-4

def fractional_flow(water_saturation):
    ff_value = krw(water_saturation)/water_viscosity
    ff_value /= ff_value+kro(water_saturation)/oil_viscosity
    return ff_value

sol = bl.buckley_leverett(2., 
                          fractional_flow = fractional_flow, 
                          residual_water = .2, 
                          residual_oil = .0, 
                          A=1., 
                          injection_rate=1., 
                          porosity=1.)


for x in range(10):
    print x, sol(x)
