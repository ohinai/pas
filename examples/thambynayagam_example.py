
import pas.thambynayagam as th
import numpy as np

case1 = th.Thambynayagam11_2()

case1.param["point_source_rate"] = 1.

pressure = case1.point_solution(.2, .2, .2, .1)

print pressure



