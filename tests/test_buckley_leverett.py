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

        case1.param["residual_n"] = .2
        case1.param["residual_w"] = .1
        case1.param["initial_sw"] = .1       

        def krw(sw):
            se = sw-case1.param["residual_w"]
            se /= 1.-case1.param["residual_n"]-case1.param["residual_w"]
            return .4*se**4

        def krn(sw):
            se = sw-case1.param["residual_w"]
            se /= 1.-case1.param["residual_n"]-case1.param["residual_w"]
            return (1.-se)**4
        
        case1.k_rw = krw
        case1.k_rn = krn
        
        case1.param["viscosity_w"] = 1.e-4
        case1.param["viscosity_n"] = 1.e-4
        
        case1.param["length"] = 5.
        case1. saturation_solution(1.)
        
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
