#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2014 Srinivasa Rao Zinka
# License: New BSD License.

from __future__ import division
import numpy as np
import filtertool as ft

N = 4
poles = np.array([1.2576, -0.1546 - 0.9218j, -0.1546 + 0.9218j])
# poles = np.array([])

F, P = ft.Chebyshev_gen(N, poles)
print F
print P

ft.plot_rational(F, P, x_min= -1.1, x_max=1.1, x_num=1000)

F = ft.w_to_s(F, coef_norm=True)
P = ft.w_to_s(P, coef_norm=True)

eps = 1.5; eps_R = 1
E, roots_E = ft.poly_E(eps, eps_R, F, P)

print E
print F
print P

ft.plot_mag(eps, eps_R, F, P, E, w_min= -6, w_max=6, w_num=500, dB=True,
             dB_limit= -70)
#
#M, Rs, Rl = ft.coupling_N(F, P, E, eps, eps_R)
#print 'M:', '\n', M.real
#print 'Rs:', Rs
#print 'Rl:', Rl
#
#ft.MN_to_Sparam(M, Rs, Rl, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)
#
#M_orig = M
#loop = True
#while(loop):
#    p=input ('Enter the element position p value:  ')
#    q=input('Enter the element position  q value:  ')
#    i=input('Enter the pivot i value:  ')
#    j=input('Enter the pivot j value: ')    
#    
#    tht=ft.annihilation(p,q,i,j,N,M)
#    R=ft.rotation_R(N,i,j,tht)
#    
#    T=np.dot(R,M)
#    M=np.dot(T,np.transpose(R))
#    print M        
#    ft.MN_to_Sparam(M, Rs, Rl, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)
#    
#    goback = input('If you want to goback to the original M, give 1')
#    if(goback):
#        M = M_orig
#        print M
#    loop = input('If you want to annihilate further, give 1')