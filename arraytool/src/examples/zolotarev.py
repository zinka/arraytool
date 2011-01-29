#! /usr/bin/env python
""" 
Zolotarev array distribution. (not finished yet)

References:
[1] McNamara, D. A., "Direct synthesis of optimum difference patterns for 
    discrete linear arrays using Zolotarev distributions", 
    IEE Proceedings H Microwaves, Antennas and Propagation, 1993, 140, 495-500 
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
    """Jacobi theta (\Theta or B as depicted in [1]) function."""
    q = qfrom(m=m)
    DM = ellipk(m)
    z = mp.pi * u / (2 * DM)
    theta = jtheta(n=4, z=z, q=q, derivative=0)
    return theta
def z_eta(u, m):
    """Jacobi eta (H) function [1]."""
    q = qfrom(m=m)
    DM = ellipk(m)
    z = mp.pi * u / (2 * DM)
    eta = jtheta(n=1, z=z, q=q, derivative=0)
    return eta
def z_zolotarev2(N, x, m):
    """Value of the Zolotarev polynomial in region II [1]."""    
    # argument (eq.(4),[1])
    M = -ellipk(m) / N
    snM = ellipfun('sn', u=M, m=m)    
    # evaluation of p (eq.(24),[1])
    p = mp.sqrt((snM ** 2 - x ** 2) / (m * snM ** 2 * (1 - x ** 2)))    
    # evaluation of s
    s = ellipf(p, m)    
    # evaluation of 'f(M,s,m)'
    fMsm = mp.log(z_theta(M + s, m) / z_theta(M - s, m))    
    # finally, the value of the Zolotarev polynomial  (eq.(21),[1])
    f = mp.cos(n * mp.pi) * mp.cosh((n + 0.5) * fMsm)        
    return f
def z_x123fromm(N, m):
    """Function to get x2, where Zolotarev polynomial is 
       maximum/minimum in region II [1]."""
    # argument
    M = -ellipk(m) / N
    # evaluation of 'x2' value for given N, k values
    snMM = ellipfun('sn', u= -M, m=m)
    snM = ellipfun('sn', u=M, m=m)
    cnM = ellipfun('cn', u=M, m=m)
    dnM = ellipfun('dn', u=M, m=m)
    znM = ellipe(M, m) - ellipe(m) * ellipf(M, m) / ellipk(m)
    x3 = snMM
    x1 = x3*mp.sqrt(1-m)/dnM
    x2 = x3 * mp.sqrt(1 - (cnM * znM) / (snM * dnM))  
    return (x1, x2, x3)
#==============================================================================
# Checking SLR value
#==============================================================================
# various parameters
n = 9
N = 2 * n + 1
k = mpf('0.999895316')
m = k ** 2 # MODULUS 'k' is not convenient, so we use PARAMETER 'm' instead
#  x1, x2 and x3 positions [1]
x1, x2, x3 = z_x123fromm(N, m)
# value of Zolotarev polynomial in region II
R = z_zolotarev2(N, x2, m)
print R
print 10 * mp.log10(R ** 2) # SLR depends only on the magnitude of R here