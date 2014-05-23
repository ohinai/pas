#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pas
----------------------------------

Tests for `buckley_leverett` module.
"""

import unittest

import pas.buckley_leverett as bl

class TestBuckleyLeverett(unittest.TestCase):

    def setUp(self):
        pass

    def test_bl(self):

        case1 = bl.BuckleyLeverett()

        krw=lambda se:se
        kro=lambda se:(1.-se)

        water_viscosity = 1.e-4
        oil_viscosity = 1.e-4

        def fractional_flow(water_saturation):
            ff_value = krw(water_saturation)/water_viscosity
            ff_value /= ff_value+kro(water_saturation)/oil_viscosity
            return ff_value

        case1.fractional_flow = fractional_flow
        sol = case1.water_saturation_solution(2.)
        
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
