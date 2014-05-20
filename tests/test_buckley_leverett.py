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
        sol = bl.buckley_leverett(1., 
                                  krw=lambda se:se, 
                                  kro=lambda se:(1.-se), 
                                  residual_water = .2, 
                                  residual_oil = .0, 
                                  water_viscosity=1.e-4, 
                                  oil_viscosity=1.e-4, 
                                  A=1., 
                                  injection_rate=1., 
                                  porosity=1.)
        
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
