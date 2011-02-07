#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

""" 
Zolotarev array distribution. (almost done!)

References:
[1] McNamara, D. A., "Direct synthesis of optimum difference patterns for 
    discrete linear arrays using Zolotarev distributions", 
    IEE Proceedings H Microwaves, Antennas and Propagation, 1993, 140, 495-500
[2] Abramowitz and Stegun, "Handbook of Mathematical Functions", Dover, 1964,
    http://people.math.sfu.ca/~cbm/aands/frameindex.htm, Page 577
[3] Levy, R., "Generalized Rational Function Approximation in Finite Intervals
    Using Zolotarev Functions", #IEEE_J_MTT#, 1970, 18, 1052-1064
"""

from __future__ import division
from mpmath import mpf, mp, qfrom, ellipk, ellipe, ellipf, ellipfun, jtheta

#==============================================================================
# Precision level
#==============================================================================

mp.dps = 25; mp.pretty = True

#==============================================================================
# Custom functions
#==============================================================================

def z_theta(u, m):
    """Jacobi theta function (eq 16.31.1, [2, p.577])."""
    q = qfrom(m=m)
    DM = ellipk(m)
    z = mp.pi * u / (2 * DM)
    theta = jtheta(n=4, z=z, q=q)
    return theta

def z_eta(u, m):
    """Jacobi eta function (eq 16.31.3, [2, p.577])."""
    q = qfrom(m=m)
    DM = ellipk(m)
    z = mp.pi * u / (2 * DM)
    eta = jtheta(n=1, z=z, q=q)
    return eta

def z_am(u, m):
    """Jacobi amplitude function (eq 16.1.5, [2, p.569])."""
    snM = ellipfun('sn', u=u, m=m)
    cnM = ellipfun('cn', u=u, m=m)
    if (0 <= cnM <= 1):
        phi = mp.asin(snM)
    elif (-1 <= cnM < 0):
        if (snM >= 0):
            phi = mp.pi - mp.asin(snM)
        else:
            phi = mp.asin(snM) - mp.pi
    else:
        print "This function only handles real 'phi' values."
    return phi

def z_zn(u, m):
    """Jacobi Zeta (zn(u,m)) function (eq 16.26.12, [2, p.576])."""
    phi = z_am(u, m)
    zn = ellipe(phi, m) - ellipe(m) * u / ellipk(m)
    return zn

def z_x123from_m(N, m):
    """Function to get x1, x2 and x3 (eq 3, 5 and 6, [1])."""
    M = -ellipk(m) / N
    snMM = ellipfun('sn', u= -M, m=m)
    snM = ellipfun('sn', u=M, m=m)
    cnM = ellipfun('cn', u=M, m=m)
    dnM = ellipfun('dn', u=M, m=m)
    znM = z_zn(M, m)
    x3 = snMM
    x1 = x3 * mp.sqrt(1 - m) / dnM
    x2 = x3 * mp.sqrt(1 - (cnM * znM) / (snM * dnM))  
    return x1, x2, x3

def z_zolotarev(N, x, m):
    """Function to evaluate the Zolotarev polynomial (eq 1, [1])."""
    M = -ellipk(m) / N
    x3 = ellipfun('sn', u= -M, m=m)  
    xbar = x3 * mp.sqrt((x ** 2 - 1) / (x ** 2 - x3 ** 2)) # rearranged eq 21, [3]
    u = ellipf(mp.asin(xbar), m) # rearranged eq 20, [3], asn(x) = F(asin(x)|m)     
    f = mp.cosh((n + 0.5) * mp.log(z_eta(M + u, m) / z_eta(M - u, m)))
    if (f.imag / f.real > 1e-10):
        print "imaginary part of the Zolotarev function is not negligible!"
        print "f_imaginary = ",f.imag
    else:
        if (x>0): # no idea why I am doing this ... anyhow, it seems working
            f = -f.real  
        else:
            f = f.real        
    return f

def z_str2num(s):
    """Convert string to either int or float."""
    try:
        ret = int(s)
    except ValueError:
        #Try float.
        ret = float(s)
    return ret

#==============================================================================
# '__main__' function
#==============================================================================

if __name__ == '__main__':        
    
    #==========================================================================
    # Basic parameters
    #==========================================================================
    
    n = 5
    N = 2 * n + 1
    k = mpf('0.999895316') # for 25 dB 
#    k = mpf('0.999971042') # for 30 dB
    m = k ** 2 # MODULUS 'k' is not convenient, so we use PARAMETER 'm' instead
    
    #==========================================================================
    # Testing side-lobe ratio (SLR)
    #==========================================================================
    
    x1, x2, x3 = z_x123from_m(N, m)
    #print x1, '\n', x2, '\n', x3

    x = x2
    R = z_zolotarev(N, x, m)
    print R
    print 10 * mp.log10(R ** 2) # SLR depends only on the magnitude of R here
    
    #==========================================================================
    # Plotting the actual polynomial
    #==========================================================================
    
    import numpy as np
    x = np.linspace(-1.04, 1.04, num=100, endpoint=True, retstep=False)
    y = []
    for i in range(len(x)):
        tmp = mp.nstr(z_zolotarev(N, x[i], m), n=6)
        tmp = z_str2num(tmp)
        y.append(tmp)
    
#    x11 = z_str2num(mp.nstr(x1, n=6))
#    x22 = z_str2num(mp.nstr(x2, n=6))
#    x33 = z_str2num(mp.nstr(x3, n=6))
    
    # Plotting
    import matplotlib.pyplot as plt
    f1 = plt.figure(1)
    p1 = plt.subplot(111)
    p1.plot(x, y)
    p1.axis('tight')
    p1.grid(True)
    plt.title(r'$Z_{2n+1}(x)$')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    p1.axhspan(-1, 1, facecolor='y', alpha=0.5)
    plt.show()