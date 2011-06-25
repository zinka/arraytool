#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

"""
Most of the times, we don't have to go through the process given in the last 
two sections. Instead, we can get the excitation values in one single step using
the function :func:`dist`. This function gives array excitation coefficients 
corresponding to various array distribution types such as Dolph-Chebyshev, 
McNamara-Zolotarev-sum, McNamara-Zolotarev-diff-f, McNamara-Zolotarev-diff-s, 
Taylor, Bayliss, Pritchard-Chebyshev-be, Pritchard-Chebyshev-ue, etc. A simple
example to obtain Taylor n-bar (also known as `Villeneuve`) distribution is shown below::

    import arraytool.planar as planar
    
    a = 0.5 # separation between the elements along x-axis (normalized WRS wavelength)
    M = 10 # no. of elements along x-axis
    
    SLR = 25 # side-lobe ratio in dB
    R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
    
    A = planar.dist(a, M, R, dist_type_x="Taylor", mbar=2, alpha_x=0)
    print 'array coefficients:', '\\n', A.T
    
The output is as shown below::    
    
    >>> arrayfactor zeros: 
    array([[ 0.52917308],
           [ 0.61909302],
           [ 0.76458654],
           [ 0.91008006],
           [ 1.        ],
           [ 1.        ],
           [ 0.91008006],
           [ 0.76458654],
           [ 0.61909302],
           [ 0.52917308]])

.. note::
           
    However, the function :func:`dist` will not provide any information regarding the
    array factor zeros. If this information is necessary, then there is no other option 
    but to use the function :func:`AF_zeros`
"""

#import arraytool.planar as planar
#
#a = 0.5 # separation between the elements along x-axis (normalized WRS wavelength)
#M = 10 # no. of elements along x-axis
#
#SLR = 25 # side-lobe ratio in dB
#R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
#
#A = planar.dist(a, M, R, dist_type_x="Taylor", mbar=2, alpha_x=0)
#print 'array coefficients:', '\n', A.T
