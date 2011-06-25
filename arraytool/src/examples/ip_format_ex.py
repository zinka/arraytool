#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

"""
One way to generate Arraytool's input format is to use the function :func:`ip_format`::

    import arraytool.planar as planar
    import numpy as np
    
    # Array lattice parameters (a & b are normalized with respect to the wavelength)
    a = 0.5 # separation between the elements along x-axis
    b = 0.7 # separation between the elements along y-axis
    gamma = np.pi / 2 # lattice angle in radians
    
    # Array Excitation information
    M = 1 # no. of elements along x-axis
    N = 11 # no. of elements along y-axis
    A = np.ones((N, M)) # A simple uniform planar excitation
    
    # Generating Arraytool input format 
    array_ip = planar.ip_format(a, b, A, gamma)
    
    print array_ip
    
In the above coding, matrix A represents the planar array arrangement which is
shown below.

TBD

Also, if needed, one can plot the array excitation using the option ``plot`` of
:func:`ip_format`. If the array is linear (either along x or y axis), a 2D plot
is generated using `Matplotlib`. If the given array is planar, a 3D stem plot is
generated using `Mayavi`. Two such simple plots are shown below.

TBD

Another way is, to simply make a CSV file in the below format and and import it 
using the function :func:`at_import`.
"""

#import arraytool.planar as planar
#import numpy as np
#
## Array lattice parameters (a & b are normalized with respect to the wavelength)
#a = 0.5 # separation between the elements along x-axis
#b = 0.7 # separation between the elements along y-axis
#gamma = np.pi / 2 # lattice angle in radians
#
## Array Excitation information
#M = 1 # no. of elements along x-axis
#N = 11 # no. of elements along y-axis
#A = np.ones((N, M)) # A simple uniform planar excitation
#
## Generating Arraytool input format 
#array_ip = planar.ip_format(a, b, A, gamma)
#
#print array_ip