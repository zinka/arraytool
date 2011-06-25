#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

"""
Array factor of a `linear discrete` array (neglecting any normalization factor)
is defined as

.. math::

  AF(u) = \sum_{m=1}^{M}\ [A_{m}\exp(jk_{0}ux_{m}].

Some of the array factor properties are given as:

- If all the elements in the array are placed uniformly, then :math:`AF\left(u\\right)`
  is a periodic function in the :math:`u` - domain with a period of :math:`\lambda/a`,
  where :math:`a` is the separation between elements.
  
- Since :math:`AF\left(u\\right)` is a polynomial of the order :math:`M-1` in :math:`u`,
  if we know :math:`M-1` zeros of the array factor (with in a period), we can easily evaluate all 
  array coefficients :math:`A_m`. Fortunately, in most of the cases (except 
  shaped beam synthesis), these zeros are symmetrically located on the :math:`u`
  axis. So, in symmetric cases, we need only :math:`\mathrm{ceil}\left[(M-1)/2\\right)]`
  zeros. If the pattern is a difference pattern, then zero at the origin also should be
  taken into consideration.
  
In Arraytool, the function :func:`AF_zeros` provides these :math:`\mathrm{ceil}\left[(M-1)/2\\right)]`
zeros. So, let us get the array factor zeros for some simple array excitations::

    import arraytool.planar as planar
    
    a = 0.5 # separation between the elements along x-axis (normalized WRS wavelength)
    M = 10 # no. of elements along x-axis
    
    SLR = 25 # side-lobe ratio in dB
    R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
    
    U0 = planar.AF_zeros(a, M, R, dist_type="Dolph-Chebyshev")
    print 'arrayfactor zeros:', '\\n', U0
    
The output is as shown below::    
    
    >>> arrayfactor zeros: 
    array([[ 0.2348126 ],
           [ 0.38767709],
           [ 0.58329742],
           [ 0.78998547]])
   
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
#print 'arrayfactor zeros:', '\n', U0