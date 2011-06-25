#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

"""
When the given array is planar, i.e., a 2D array in the :math:`xy`-plane, then 
we need to plot the pattern as a function of :math:`(u,v)`. In order to accomplish
this task, one can use the function :func:`pattern_uv`. A simple example to plot
radiation pattern of a planar uniform array is shown below::

    import arraytool.planar as planar
    import numpy as np
    
    a = 0.6 # separation between the elements along x-axis (normalized WRS wavelength)
    b = 0.5 # separation between the elements along y-axis (normalized WRS wavelength)
    M = 10 # no. of elements along x-axis
    N = 11 # no. of elements along y-axis
    
    A = np.ones((N, M)) # Uniform planar excitation
    
    # Converting the 'excitation & position' information into Arraytool input format
    array_ip = planar.ip_format(a, b, A)
    
    # Calling the 'pattern_uv' function to evaluate and plot 3D AF/GF/NF
    planar.pattern_uv(array_ip, u_scan=0, v_scan=0, u_min= -2, u_max=2, u_num=300,
               v_min= -2, v_max=2, v_num=300, scale="dB", dB_limit=-40,
               factor="NF", plot_type="rect", mayavi_app=False)
               
If a contour plot is needed, then the following command can be used::

    # Calling the 'pattern_uv' function to evaluate and plot 3D AF/GF/NF
    planar.pattern_uv(array_ip, u_scan=0, v_scan=0, u_min= -2, u_max=2, u_num=300,
               v_min= -2, v_max=2, v_num=300, scale="dB", dB_limit=-40,
               factor="NF", plot_type="contour", mayavi_app=False)              
"""

#import arraytool.planar as planar
#import numpy as np
#
#a = 0.6 # separation between the elements along x-axis (normalized WRS wavelength)
#b = 0.5 # separation between the elements along y-axis (normalized WRS wavelength)
#M = 10 # no. of elements along x-axis
#N = 11 # no. of elements along y-axis
#
#A = np.ones((N, M)) # Uniform planar excitation
#
## Converting the 'excitation & position' information into Arraytool input format
#array_ip = planar.ip_format(a, b, A)
#
## Calling the 'pattern_uv' function to evaluate and plot 3D AF/GF/NF
#planar.pattern_uv(array_ip, u_scan=0, v_scan=0, u_min= -2, u_max=2, u_num=300,
#           v_min= -2, v_max=2, v_num=300, scale="dB", dB_limit=-40,
#           factor="NF", plot_type="rect", mayavi_app=False)