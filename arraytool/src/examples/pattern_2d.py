#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

""" Module for 2D PATTERN EVALUATION .... """

from __future__ import division
from numpy import ones, mgrid, array, zeros_like, hstack, pi, \
                  reshape, tan, random, angle, flipud, sqrt, set_printoptions, \
                  nan, empty, sin 
from enthought.mayavi import mlab
from enthought.tvtk.api import tvtk

set_printoptions(precision=2, threshold=nan, edgeitems=None,
                 linewidth=None, suppress=None, nanstr=None, infstr=None)

#==============================================================================
# Function "ip_format"
#==============================================================================

def ip_format(a, b, Amn, gamma=pi / 2):
    """ Function to generate array input data (in specified format) """
    
    M = float(Amn.shape[1])    # don't get confused here ... (numrows, numcols = X.shape)
    N = float(Amn.shape[0])    # M, N = No. of columns & rows of physical array

    if (M == 1):  # i.e, linear array along y-direction
        a = 0
        gamma = pi / 2
        
    if (N == 1):  # i.e, linear array along x-direction
        b = 0
        gamma = pi / 2    
        
    xlim = (M * a) / 2 
    ylim = ((N - 1) * b) / 2  
        
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
            
    x = reshape(x, (M * N, -1))
    y = reshape(y, (M * N, -1))
    z = zeros_like(x)   # because only planar arrays on XY-plane are considered
    Amn = reshape(Amn, (M * N, -1))
    array_ip = hstack((x, y, z, Amn))  # array excitation IMPORT format
    
    return array_ip

#==============================================================================
# Function "pattern_2d"
#==============================================================================

def pattern_2d(array_ip, freq, scan_freq, tht_resol, phi_resol,
               scan_tht=0, scan_phi=0, m_bits=0, n_bits=0,
               tht_cut=pi / 2, phi_cut=0,
               tht_min= -pi, tht_max=pi, phi_min= -pi, phi_max=pi,
               FFT=False, planar=True):
    """" Function to evaluate 2D-pattern cuts (Theta & Phi cuts) """
    
    x = array_ip[:, 0]
    y = array_ip[:, 1]
    z = array_ip[:, 2]
    Amn = array_ip[:, 3]
    dis = sqrt(x ** 2 + y ** 2)   # do the "pick" operation later    
    
    #==========================================================================
    # Amplitude and phase quantizations
    #==========================================================================    
    
    norm_Amn = abs(Amn).max()
    Amn = Amn / norm_Amn
    
    B = abs(Amn)
    PH = angle(Amn)
    
    #==========================================================================
    # Plotting 3D-plots
    #==========================================================================
#    mlab.options.backend = 'envisage'         # one way to save visualization
#    f = mlab.figure()
    
#    s1 = mlab.quiver3d(x, y, z, z, z, Amn, scale_factor=1)   # stem representation  
    s1 = mlab.quiver3d(x, y, z, z, z, Amn)   # stem representation  
#    mlab.points3d(x, y, z, scale_factor=0.005)

    #==========================================================================
    # defining outline "ranges" and plotting the XYZ-coordinate system
    #==========================================================================
    
    xmin = min(x)
    print xmin
    xmax = max(x)
    ymin = min(y)
    ymax = max(y)
    zmin = min(Amn)
    zmax = max(Amn)
    ranges1 = [xmin, xmax, ymin, ymax, zmin, zmax]
    ranges1 = [min(x), max(x), min(y), max(y), min(Amn), max(Amn)]
    
#    mlab.outline()
    a1 = mlab.axes(xlabel="X", ylabel="Y", zlabel="Amn", ranges=ranges1)
    zmax = a1.axes.bounds[5]
    
    #==========================================================================
    # GRID generation (better to write a simple function for this)
    #==========================================================================
    
    # Generate some points.
    x, y, z = mgrid[xmin:xmax:11j, ymin:ymax:13j, 0:zmax:6j]
    
    # The actual points.
    pts = empty(z.shape + (3,), dtype=float)
    pts[..., 0] = x
    pts[..., 1] = y
    pts[..., 2] = z
    
    # We reorder the points, scalars and vectors so this is as per VTK's
    # requirement of x first, y next and z last.
    pts = pts.transpose(2, 1, 0, 3).copy()
    pts.shape = pts.size / 3, 3
    
    # Create the data-set.
    sg = tvtk.StructuredGrid(dimensions=x.shape, points=pts)
    
    # Now visualize the data.
    d = mlab.pipeline.add_dataset(sg)
    gx = mlab.pipeline.grid_plane(d)
    gy = mlab.pipeline.grid_plane(d)
    gy.grid_plane.axis = 'y'
    gz = mlab.pipeline.grid_plane(d)
    gz.grid_plane.axis = 'z'
    
    #==========================================================================
    # 
    #==========================================================================

    mlab.show_pipeline()
    s1.scene.isometric_view()
    mlab.show()    
    
#==============================================================================
# Testing all functions
#==============================================================================

M = 5     # no. of elements along x-axis
N = 6     # no. of elements along y-axis
a = 0.015    # separation along x-axis
b = 0.010    # separation along y-axis
gamma = pi / 2 # lattice angle
#Amn = ones((N, M)) # Amn should be given in complex(a+bj) form
Amn = random.rand(N, M) 
#Amn = array([[1 + 3j, 0, 1], [1, 2, 3]])

#-----------------------------------------------------------------------------
 
freq = 10e9
scan_freq = 10e9
scan_tht = (0) * (pi / 180)
scan_phi = (0) * (pi / 180)
m_bits = 0
n_bits = 0
tht_cut = (90) * (pi / 180)
phi_cut = (0) * (pi / 180)
tht_resol = (1) * (pi / 180)
phi_resol = (1) * (pi / 180)
tht_min = (-180) * (pi / 180)
tht_max = (180) * (pi / 180)
phi_min = (-180) * (pi / 180)
phi_max = (180) * (pi / 180)

#----------------------------------------------------------------------------- 

if __name__ == "__main__":
    array_ip = ip_format(a, b, Amn, gamma=pi / 2.5)
#    print array_ip
    pattern_2d(array_ip, freq, scan_freq, tht_resol, phi_resol,
               scan_tht=0, scan_phi=0, m_bits=0, n_bits=0,
               tht_cut=pi / 2, phi_cut=0,
               tht_min= -pi, tht_max=pi, phi_min= -pi, phi_max=pi, FFT=False)
    
#----------------------------------------------------------------------------- 

