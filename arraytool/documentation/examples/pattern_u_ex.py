#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

"""
Obviously, the most important task in array analysis and design is to plot the actual
array factor. In order to accomplish this task, one can use the function :func:`pattern_u`.
A simple example to plot a Dolph-Chebyshev pattern is shown below::

    import arraytool.planar as planar
    
    a = 0.6 # separation between the elements along x-axis (normalized WRS wavelength)
    M = 10 # no. of elements along x-axis
    
    SLR = 25 # side-lobe ratio in dB
    R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
    
    A = planar.dist(a, M, R, dist_type_x="Dolph-Chebyshev")
    
    # Converting the 'excitation & position' information into Arraytool input format
    array_ip = planar.ip_format(a, b=0, A=A)
    
    # Calling the 'pattern_u' function to evaluate and plot 2D AF/GF/NF
    [u,AF] = planar.pattern_u(array_ip, u_scan=0, u_min= -2, u_max=2, u_num=300,
                              scale="dB", dB_limit= -40, factor="NF", 
                              plot_type="rect", lattice=True)
                              
.. image:: _static/pattern_u_1.png                              
                              
If a polar plot is needed, then the following command can be used::

    # Calling the 'pattern_u' function to evaluate and plot 2D AF/GF/NF
    [u,AF] = planar.pattern_u(array_ip, u_scan=0, u_min= -1, u_max=1, u_num=300,
                              scale="dB", dB_limit= -40, factor="NF", 
                              plot_type="polar", lattice=True)
                              
.. image:: _static/pattern_u_2.png
"""

#import arraytool.planar as planar
#
#a = 0.6 # separation between the elements along x-axis (normalized WRS wavelength)
#M = 10 # no. of elements along x-axis
#
#SLR = 25 # side-lobe ratio in dB
#R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
#
#A = planar.dist(a, M, R, dist_type_x="Dolph-Chebyshev")
#
## Converting the 'excitation & position' information into Arraytool input format
#array_ip = planar.ip_format(a, b=0, A=A)
#
## Calling the 'pattern_u' function to evaluate and plot 2D AF/GF/NF
#[u,AF] = planar.pattern_u(array_ip, u_scan=0, u_min= -2, u_max=2, u_num=300,
#                          scale="dB", dB_limit= -40, factor="NF", 
#                          plot_type="rect", lattice=True)
#
### Calling the 'pattern_u' function to evaluate and plot 2D AF/GF/NF
##[u,AF] = planar.pattern_u(array_ip, u_scan=0, u_min= -1, u_max=1, u_num=300,
##                          scale="dB", dB_limit= -40, factor="NF", 
##                          plot_type="polar", lattice=True)