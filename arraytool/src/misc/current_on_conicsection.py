#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2014 Srinivasa Rao Zinka
# License: New BSD License.

""" A simple program to plot induced surface current on paraboloid, spheroid 
etc ... """

from __future__ import division
from numpy import sqrt, arange, pi, cos, sin, mgrid, ones_like, arccos, \
set_printoptions, nan, array, reshape, cross, hstack
#from enthought.mayavi import mlab
from mayavi import mlab

set_printoptions(precision=2, threshold=nan, edgeitems=None,
                 linewidth=None, suppress=None, nanstr=None, infstr=None)

#==============================================================================
# General ***boid 
#==============================================================================

# Rectangular data
l1 = 10
eps1 = 1
[x, y] = mgrid[-10:10:0.5, -10:10:0.5]

if (eps1 == 1):
    
    z2 = (l1 ** 2 - x ** 2 - y ** 2) / (2 * l1)
    r2 = sqrt(x ** 2 + y ** 2 + z2 ** 2)
    tht2 = arccos(z2 / r2)
    phi2 = arccos(x / (r2 * sin(tht2))) 
        
    u2 = -(x) / (sqrt(l1 ** 2 - (1 - eps1 ** 2) * (x ** 2 + y ** 2)))
    v2 = -(y) / (sqrt(l1 ** 2 - (1 - eps1 ** 2) * (x ** 2 + y ** 2)))
    w2 = -ones_like(u2) 
        
else:
    
    z1 = -(l1 * eps1) / (1 - eps1 ** 2)     \
                    - (l1) / (1 - eps1 ** 2) * sqrt(1 - ((1 - eps1 ** 2) * (x ** 2 + y ** 2)) / (l1 ** 2))
    u1 = -(x) / (sqrt(l1 ** 2 - (1 - eps1 ** 2) * (x ** 2 + y ** 2)))
    v1 = -(y) / (sqrt(l1 ** 2 - (1 - eps1 ** 2) * (x ** 2 + y ** 2)))
    w1 = +ones_like(u1)                     
                    
    z2 = -(l1 * eps1) / (1 - eps1 ** 2)     \
                    + (l1) / (1 - eps1 ** 2) * sqrt(1 - ((1 - eps1 ** 2) * (x ** 2 + y ** 2)) / (l1 ** 2))
    r2 = sqrt(x ** 2 + y ** 2 + z2 ** 2)
    tht2 = arccos(z2 / r2)
    phi2 = arccos(x / (r2 * sin(tht2)))                    
                    
    u2 = -(x) / (sqrt(l1 ** 2 - (1 - eps1 ** 2) * (x ** 2 + y ** 2)))
    v2 = -(y) / (sqrt(l1 ** 2 - (1 - eps1 ** 2) * (x ** 2 + y ** 2)))
    w2 = -ones_like(u2) 
                    
## Polar data
#[tht,phi] = mgrid[0:pi/2:pi/20,0:2*pi+pi/20:pi/20]
#
#r = l1/(1+eps1*cos(tht))
#x_r = r*sin(tht)*cos(phi)
#y_r = r*sin(tht)*sin(phi)
#z_r = r*cos(tht)

#==============================================================================
# H-field
#==============================================================================

#H_tht = sin(phi2)*(cos(tht2)+eps1)
#H_phi = cos(phi2)*(eps1*cos(tht2)+1)    # Koffman field (only minus difference)

H_tht = -sin(phi2) * (cos(tht2) + eps1)
H_phi = -cos(phi2) * (eps1 * cos(tht2) + 1)    # actual field (right now, equal to Koffman's)

H_x = cos(tht2) * cos(phi2) * H_tht - sin(phi2) * H_phi
H_y = cos(tht2) * sin(phi2) * H_tht + cos(phi2) * H_phi
H_z = -sin(tht2) * H_tht

#==============================================================================
# CROSS product evaluation
#==============================================================================

H_x_temp = reshape(H_x, (H_x.size, 1))
H_y_temp = reshape(H_y, (H_x.size, 1))
H_z_temp = reshape(H_z, (H_x.size, 1))
H_tot_temp = hstack((H_x_temp, H_y_temp, H_z_temp))

u2_temp = reshape(u2, (u2.size, 1))
v2_temp = reshape(v2, (v2.size, 1))
w2_temp = reshape(w2, (w2.size, 1))
uvw_tot_temp = hstack((u2_temp, v2_temp, w2_temp))

J_temp = cross(uvw_tot_temp, H_tot_temp)
J_x_temp = J_temp[:, 0]
J_y_temp = J_temp[:, 1]
J_z_temp = J_temp[:, 2]

J_x = reshape(J_x_temp, u2.shape)
J_y = reshape(J_y_temp, u2.shape)
J_z = reshape(J_z_temp, u2.shape)

#==============================================================================
# Plotting current vectors
#==============================================================================

#mlab.options.backend = 'envisage'         # one way to save visualization
#f = mlab.figure()
c1 = mlab.points3d(0, 0, 0)
mlab.quiver3d(x, y, z2, J_x, J_y, J_z, colormap="jet", scale_factor=1)

#==============================================================================
# Plotting actual ***boids
#==============================================================================

if (eps1 == 1):
    
    s2 = mlab.mesh(x, y, z2, colormap="gray", transparent=True)
#    mlab.quiver3d(x, y, z2, u2, v2, w2, scale_factor=1)
else:
    
#    sr = mlab.mesh(x_r, y_r, z_r, colormap="gray", transparent=True)

    s2 = mlab.mesh(x, y, z2, colormap="gray", transparent=True)
#    mlab.quiver3d(x, y, z2, u2, v2, w2,colormap="gray",transparent=True, scale_factor=1)
    
#    s1 = mlab.mesh(x, y, z1,colormap="gray", transparent=True) 
#    mlab.quiver3d(x, y, z1, u1, v1, w1,colormap="gray",transparent=True, scale_factor=1) 

#==============================================================================
# Show the scene
#==============================================================================

mlab.show()

##==============================================================================
## reshaping
##==============================================================================
#
#a = array([[1,2],[3,4]])
#
#print reshape(a, (a.size,1))
#print reshape(a, a.shape)
#
##==============================================================================
## cross product
##==============================================================================
#
#x = array([[1,2,3],[1,2,3]])
#y = array([[4,5,6],[4,5,6]])
#print cross(x,y)

