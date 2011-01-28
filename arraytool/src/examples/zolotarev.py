#! /usr/bin/env python
""" 
Zolotarev array distribution. (not finished yet)

References:
[1] McNamara, D. A., "Direct synthesis of optimum difference patterns for 
discrete linear arrays using Zolotarev distributions", 
IEE Proceedings H Microwaves, Antennas and Propagation, 1993, 140, 495-500 
"""

from __future__ import division
from mpmath import mpf, mp, qfrom, ellipk, ellipe, ellipf, ellipfun, jtheta, mpc

#==============================================================================
# Precision level
#==============================================================================
mp.dps = 25; mp.pretty = True
#==============================================================================
# Custom functions
#==============================================================================
def z_theta(u, m):
    """Jacobi theta (\Theta) function."""
    q = qfrom(m=m)
    DM = ellipk(m)
    z = mp.pi * u / (2 * DM)
    theta = jtheta(n=4, z=z, q=q, derivative=0)
    return theta
def z_eta(u, m):
    """Jacobi eta (H) function."""
    q = qfrom(m=m)
    DM = ellipk(m)
    z = mp.pi * u / (2 * DM)
    eta = jtheta(n=1, z=z, q=q, derivative=0)
    return eta
#==============================================================================
# Checking SLR value
#==============================================================================
n = 9
N = 2 * n + 1

# various parameters
k1 = mpf('0.999971042')
m1 = k1 ** 2
q1 = qfrom(m=m1)

# argument
M1 = -ellipk(m1)/N

# evaluation of 'x2' value for given N, k values
snMM = ellipfun('sn', u=-M1, m=m1)
snM = ellipfun('sn', u=M1, m=m1)
cnM = ellipfun('cn', u=M1, m=m1)
dnM = ellipfun('dn', u=M1, m=m1)
znM = ellipe(M1,m1)-ellipe(m1)*ellipf(M1,m1)/ellipk(m1)
x2 = snMM*mp.sqrt(1-(cnM*znM)/(snM*dnM))

# evaluation of p2
p2 = mp.sqrt((snM**2-x2**2)/(m1*snM**2*(mpf('1')-x2**2)))

# evaluation of s2
s2 = ellipf(p2, m1)

# evaluation of 'f'
f = mp.log(z_theta(M1+s2, m1)/z_theta(M1-s2, m1))

# evaluation of R
R = mp.cos(n*mp.pi)*mp.cosh((n+0.5)*f)
#R = mp.cosh((n+0.5)*f)

print R
print 10*mp.log10(R**2) # SLR depends only on the magnitude of R here
#==============================================================================
# checking Jacobi THETA (\Theta) and ETA (H) functions
#==============================================================================
#M1 = 0.7
#m1 = 0.5
#theta = z_theta(M1,m1)
#eta = z_eta(M1,m1)
#print theta, '\n', eta