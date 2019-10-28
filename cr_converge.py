# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:33:59 2019

@author: Jahmour
"""
#Section of code designed to make theory and model correlation functions agree 

import numpy

N = 21
avals = [alphaV,alphaH,alphaD]
avals_nl = [0,0,0]
sigma_a = 1
if True:
    tol = 1.e-5 #Pick a tolerance below which the two Crs are considered equal
    fsBFE_out = 2*sBFE_out+1
    BFEK_model = numpy.zeros((fsBFE_out,fsBFE_out))
    element_diff = 10
    while element_diff > tol:
        theory_Cr = solve_corr(BFEK_model,N,I,gain,beta,sigma_a,tslices,avals,avals_nl)\
        *((g**2)/(I**2*(tb-ta)*(td-tc)))
        observed_Cr = BFEK
        difference = theory_Cr - observed_Cr
        element_diff = numpy.amax(abs(difference))
        BFEK_model -= difference