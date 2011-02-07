#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

""" Module for PLANAR ARRAY ANTENNAS ...."""

from __future__ import division
from numpy import mgrid, zeros_like, hstack, pi, reshape, tan, flipud
#from enthought.mayavi import mlab

#==============================================================================
# Function "ip_format"
#==============================================================================

def ip_format(a, b, Amn, gamma=pi / 2):
    """ Function to generate array input data (in specified format) """
    
    M = float(Amn.shape[1])    # don't get confused here ... 
    N = float(Amn.shape[0])    # M = No. of columns & N = No. of rows (of physical array)
    
    if (M == 1):  # i.e, linear array along y-direction
        a = 0
        gamma = pi / 2
        
    if (N == 1):  # i.e, linear array along x-direction
        b = 0
        gamma = pi / 2    
        
    xlim = (M * a) / 2 
    ylim = (N * b) / 2  
        
    #==========================================================================
    # Grid generation
    #==========================================================================
    
    [x, y] = mgrid[0:M, 0:N]
    x = (x - (M - 1) / 2).T
    y = flipud((y - (N - 1) / 2).T)
    x = x * a
    y = y * b
    
    if (gamma != pi / 2):
        x = x + (y / tan(gamma)) % a
        
    for i1 in range(int(N)):
        
        if (x[i1, 0] < -xlim):
            x[i1, :] = x[i1, :] + a
            
        if (x[i1, -1] > xlim):
            x[i1, :] = x[i1, :] - a   
            
    #==========================================================================
    # Array input data generation
    #==========================================================================
            
    x = reshape(x, (M * N, -1))
    y = reshape(y, (M * N, -1))
    z = zeros_like(x)  # because only planar arrays on XY-plane are considered
    Amn = reshape(Amn, (M * N, -1))
    array_ip = hstack((x, y, z, Amn))  # array excitation IMPORT format
    return array_ip

