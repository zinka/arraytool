#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

"""
In the previous section we obtained the array factor zeros. Now, we can use those 
zeros to obtain the array excitation using the function :func:`A_frm_zeros` as 
shown below::

    import arraytool.planar as planar
    
    a = 0.5 # separation between the elements along x-axis (normalized WRS wavelength)
    M = 10 # no. of elements along x-axis
    
    SLR = 25 # side-lobe ratio in dB
    R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
    
    U0 = planar.AF_zeros(a, M, R, dist_type="Dolph-Chebyshev")
    
    A = planar.A_frm_zeros(U0, a, M, symmetry="even").T # finding excitation coefficients
    print 'array coefficients:', '\n', A.T
    
The output is as shown below::    
    
    >>> array coefficients:
    array([[ 0.64163439],
           [ 0.59442917],
           [ 0.77799478],
           [ 0.921367  ],
           [ 1.        ],
           [ 1.        ],
           [ 0.921367  ],
           [ 0.77799478],
           [ 0.59442917],
           [ 0.64163439]])

As can be seen above, array coefficients at the center are normalized to the value 1.
"""

#import arraytool.planar as planar
#
#a = 0.5 # separation between the elements along x-axis (normalized WRS wavelength)
#M = 10 # no. of elements along x-axis
#
#SLR = 25 # side-lobe ratio in dB
#R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
#
#U0 = planar.AF_zeros(a, M, R, dist_type="Dolph-Chebyshev")
#
#A = planar.A_frm_zeros(U0, a, M, symmetry="even").T # finding excitation coefficients
#print 'array coefficients:', '\\n', A.T