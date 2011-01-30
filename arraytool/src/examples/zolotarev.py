#! /usr/bin/env python

""" 
Zolotarev array distribution. (not finished yet)

References:
[1] McNamara, D. A., "Direct synthesis of optimum difference patterns for 
    discrete linear arrays using Zolotarev distributions", 
    IEE Proceedings H Microwaves, Antennas and Propagation, 1993, 140, 495-500 
"""

from __future__ import division
from mpmath import mpf, mp, qfrom, ellipk, ellipe, ellipf, ellipfun, jtheta, \
                   nsum, inf, elliprf

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
    phi = mp.asin(snM) # decide what arcsin u need to take ... understand Jacobi amplitude
    znM = ellipe(phi, m) - ellipe(m) * ellipf(phi, m) / ellipk(m)
    x3 = snMM
    x1 = x3 * mp.sqrt(1 - m) / dnM
    x2 = x3 * mp.sqrt(1 - (cnM * znM) / (snM * dnM))  
    return x1, x2, x3

def z_f(a, b, m):
    """Summation of the series (eq.(23),[1])"""
    K = ellipk(m)
    K1 = ellipk(1 - m)
    q1 = mp.exp(-mp.pi * K / K1)
    NM = mp.cosh((a + b) * mp.pi / (2 * K))
    DM = mp.cosh((a - b) * mp.pi / (2 * K))
    fun = lambda r: ((-1) ** r / r) * ((q1 ** (2 * r)) / (1 - q1 ** (2 * r)))\
     * mp.sinh(mp.pi * r * a / K1) ** mp.sinh(mp.pi * r * b / K1)
    f = -(mp.pi * a * b / (K * K1)) + mp.log(NM / DM) - 4 * nsum(fun, [1, inf])
    return f

def z_zolotarev2(N, x, m):
    """Value of the Zolotarev polynomial in region II [1]."""    
    # argument (eq.(4),[1])
    M = -ellipk(m) / N
    snM = ellipfun('sn', u=M, m=m)    
    # evaluation of p (eq.(24),[1])
    p = mp.sqrt((snM ** 2 - x ** 2) / (m * snM ** 2 * (1 - x ** 2))) 
    # evaluation of s
    s = ellipf(p, m) # convert p to Jacobi amplitude 'phi'  and # decide what arcsin u need to take ... understand Jacobi amplitude
    # evaluation of 'f(M,s,m)'
    fMsm = mp.log(z_theta(M + s, m) / z_theta(M - s, m))    
    # finally, the value of the Zolotarev polynomial  (eq.(21),[1])
    f = mp.cos(n * mp.pi) * mp.cosh((n + 0.5) * fMsm)        
    return f

#==============================================================================
# Checking SLR value
#==============================================================================

# various parameters
n = 9
N = 2 * n + 1
k = mpf('0.999895316')
#k = mpf('0.9')
m = k ** 2 # MODULUS 'k' is not convenient, so we use PARAMETER 'm' instead
#  x1, x2 and x3 positions [1]
x1, x2, x3 = z_x123fromm(N, m)
x = x2
# value of Zolotarev polynomial in region II
R = z_zolotarev2(N, x, m)
print R
print 10 * mp.log10(R ** 2) # SLR depends only on the magnitude of R here

#==============================================================================
# Checking function f
#==============================================================================

# argument (eq.(4),[1])
M = -ellipk(m) / N
#print M
## ellipk in-terms of elliprf
#M = -elliprf(0, 1-m, 1) / N
#print M
snM = ellipfun('sn', u=M, m=m)    
# evaluation of p (eq.(24),[1])
p = mp.sqrt((snM ** 2 - x ** 2) / (m * snM ** 2 * (1 - x ** 2)))    
p = mp.acos(p)
# evaluation of s
p = 4
m = mpf(2)/3
snp = ellipfun('sn', u=p, m=m)   
phi_p = 1*mp.pi - mp.asin(snp) # decide what arcsin u need to take
print phi_p
s = ellipf(p, m)  

fMsm = z_f(M, s, m)
R = mp.cos(n * mp.pi) * mp.cosh((n + 0.5) * fMsm)        
print R
#print 10 * mp.log10(R ** 2) # SLR depends only on the magnitude of R here