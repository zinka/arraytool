#! /usr/bin/env python

""" 
Zolotarev array distribution. (not finished yet)

References:
[1] McNamara, D. A., "Direct synthesis of optimum difference patterns for 
    discrete linear arrays using Zolotarev distributions", 
    IEE Proceedings H Microwaves, Antennas and Propagation, 1993, 140, 495-500
[2] Abramowitz and Stegun, "Handbook of Mathematical Functions", Dover, 1964,
    http://people.math.sfu.ca/~cbm/aands/frameindex.htm, Page 577
"""

from __future__ import division
from mpmath import mpf, mp, qfrom, ellipk, ellipe, ellipf, ellipfun, jtheta, \
                   nsum, inf, mpc

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
    if (0<=cnM<=1):
        phi = mp.asin(snM)
    elif (-1<=cnM<0):
        if (snM>=0):
            phi = mp.pi - mp.asin(snM)
        else:
            phi = mp.asin(snM) - mp.pi
    else:
        print "This function only handles real 'phi' values."
    return phi

def z_zn(u, m):
    """Jacobi Zeta (zn(u,m)) function (eq 16.26.12, [2, p.576])."""
    phi = z_am(u, m)
    zn = ellipe(phi, m) - ellipe(m)*u/ellipk(m)
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

def z_zolotarev2(N, x, m):
    """Value of the Zolotarev polynomial in region II [1]."""    
    # argument (eq.(4),[1])
    M = -ellipk(m) / N
    snM = ellipfun('sn', u=M, m=m)  
    # evaluation of p (eq 24,[1])
    p = mp.sqrt((snM ** 2 - x ** 2) / (m * snM ** 2 * (1 - x ** 2))) 
    p = mp.asin(p) # since the notation used in [1] is a little bit different
    # evaluation of s
    s = ellipf(p, m)
    # evaluation of 'f(M,s,m)'
    fMsm = mp.log(z_theta(M + s, m) / z_theta(M - s, m))    
    # finally, the value of the Zolotarev polynomial  (eq 21,[1])
    f = mp.cos(n * mp.pi) * mp.cosh((n + 0.5) * fMsm)        
    return f

#==============================================================================
# Basic parameters
#==============================================================================

n = 9
N = 2 * n + 1
k = mpf('0.999895316')
#k = mpf('0.9')
m = k ** 2 # MODULUS 'k' is not convenient, so we use PARAMETER 'm' instead

#==============================================================================
# Checking SLR value
#==============================================================================

#  x1, x2 and x3 positions [1]
x1, x2, x3 = z_x123from_m(N, m)
print x1, '\n', x2, '\n', x3
x = x2
# value of Zolotarev polynomial in region II
R = z_zolotarev2(N, x, m)
print R
print 10 * mp.log10(R ** 2) # SLR depends only on the magnitude of R here

#==============================================================================
# Checking function f
#==============================================================================

## argument (eq.(4),[1])
#M = -ellipk(m) / N
##print M
### ellipk in-terms of elliprf
##M = -elliprf(0, 1-m, 1) / N
##print M
#snM = ellipfun('sn', u=M, m=m)    
## evaluation of p (eq.(24),[1])
#p = mp.sqrt((snM ** 2 - x ** 2) / (m * snM ** 2 * (1 - x ** 2)))    
#p = mp.acos(p)
## evaluation of s
#p = 4
#m = mpf(2)/3
#snp = ellipfun('sn', u=p, m=m)   
#phi_p = 1*mp.pi - mp.asin(snp) # decide what arcsin u need to take
#print phi_p
#s = ellipf(p, m)  
#
#fMsm = z_f(M, s, m)
#R = mp.cos(n * mp.pi) * mp.cosh((n + 0.5) * fMsm)        
#print R
##print 10 * mp.log10(R ** 2) # SLR depends only on the magnitude of R here