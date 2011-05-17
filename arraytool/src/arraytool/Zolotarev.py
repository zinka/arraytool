#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

""" 
Zolotarev polynomial related routines.

**References**

- [McNamara93]_
- [Abramowitz]_
- [Levy70]_

.. For the full description of the above citations, refer to the file references.py
    
"""

from __future__ import division
import numpy as np
from scipy import optimize
from mpmath import mpf, mp, qfrom, ellipk, ellipe, ellipf, ellipfun, jtheta
import matplotlib as mpl

# 'mpmath' precision level
mp.dps = 25; mp.pretty = True

# changing some of the 'Matplotlib' parameters
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'

def z_str2num(s):
    """Convert a string to either int or float."""
    try:
        ret = int(s)
    except ValueError:
        #Try float.
        ret = float(s)
    return ret

def z_theta(u, m):
    """Jacobi theta function (eq 16.31.1, [Abramowitz]_)."""
    q = qfrom(m=m)
    DM = ellipk(m)
    z = mp.pi * u / (2 * DM)
    theta = jtheta(n=4, z=z, q=q)
    return theta

def z_eta(u, m):
    """Jacobi eta function (eq 16.31.3, [Abramowitz]_)."""
    q = qfrom(m=m)
    DM = ellipk(m)
    z = mp.pi * u / (2 * DM)
    eta = jtheta(n=1, z=z, q=q)
    return eta

def z_am(u, m):
    """Jacobi amplitude function (eq 16.1.5, [Abramowitz]_)."""
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
    """Jacobi Zeta (zn(u,m)) function (eq 16.26.12, [Abramowitz]_)."""
    phi = z_am(u, m)
    zn = ellipe(phi, m) - ellipe(m) * u / ellipk(m)
    return zn

def z_x123_frm_m(N, m):
    """Function to get x1, x2 and x3 (eq 3, 5 and 6, [McNamara93]_)."""
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

def z_Zolotarev(N, x, m):
    """Function to evaluate the Zolotarev polynomial (eq 1, [McNamara93]_)."""
    M = -ellipk(m) / N
    x3 = ellipfun('sn', u= -M, m=m)  
    xbar = x3 * mp.sqrt((x ** 2 - 1) / (x ** 2 - x3 ** 2)) # rearranged eq 21, [Levy70]_
    u = ellipf(mp.asin(xbar), m) # rearranged eq 20, [Levy70]_, asn(x) = F(asin(x)|m)     
    f = mp.cosh((N / 2) * mp.log(z_eta(M + u, m) / z_eta(M - u, m)))
    if (f.imag / f.real > 1e-10):
        print "imaginary part of the Zolotarev function is not negligible!"
        print "f_imaginary = ", f.imag
    else:
        if (x > 0): # no idea why I am doing this ... anyhow, it seems working
            f = -f.real  
        else:
            f = f.real        
    return f

def z_Zolotarev_x2(N, m):
    """
    This function evaluates the Zolotarev polynomial at x2 for given 'N' and 'm'.
    """    
    x = z_x123_frm_m(N, m)
    ret = z_Zolotarev(N, x[1], m)
    return ret

def z_m_frm_R(N, R, a=0.1, b=0.9999999999999):
    """
    Function to obtain the parameter 'm' for a given 'SLR' and order of the 
    polynomial 'N'.
    """
    fun = lambda m: abs(z_Zolotarev_x2(N, m)) - abs(R)
    m = optimize.brentq(fun, a, b)
    return m

def z_Zolotarev_poly(N, m, interp_num=100, full=False):
    """
    Function to evaluate the polynomial coefficients of  the Zolotarev
    polynomial.
    """
    m = mpf(m)
    # Evaluating the polynomial at some discrete points    
    x = np.linspace(-1.04, 1.04, interp_num)
    y = []
    for i in range(len(x)):
        tmp = mp.nstr(z_Zolotarev(N, x[i], m), n=6)
        tmp = z_str2num(tmp)
        y.append(tmp)
    # Interpolating the data to fit to a Nth order polynomial
    coef = np.polyfit(x, y, N, full)
    roots = np.sort(np.roots(coef))
    return coef, roots

def Zolotarev2(p, q, m):
    """Another way to generate Zolotarev polynomial. It is not done yet."""    
    n = p + q
    u0 = (p / (p + q)) * ellipk(m)
    wp = 2 * (ellipfun('cd', u0, m)) ** 2 - 1
    ws = 2 * (ellipfun('cn', u0, m)) ** 2 - 1
    wq = (wp + ws) / 2
    sn = ellipfun('sn', u0, m)
    cn = ellipfun('cn', u0, m)
    dn = ellipfun('dn', u0, m)    
    wm = ws + 2 * z_zn(u0, m) * ((sn * cn) / (dn))
    beta = np.zeros((n + 5, 1))
    beta[n] = 1    
    d = np.zeros((6, 1))
    # Implementation of the main recursive algorithm
    for m1 in range(n + 2, 2, -1):        
        d[0] = (m1 + 2) * (m1 + 1) * wp * ws * wm
        d[1] = -(m1 + 1) * (m1 - 1) * wp * ws - (m1 + 1) * (2 * m1 + 1) * wm * wq
        d[2] = wm * (n ** 2 * wm ** 2 - m1 ** 2 * wp * ws) + m1 ** 2 * (wm - wq) + 3 * m1 * (m1 - 1) * wq
        d[3] = (m1 - 1) * (m1 - 2) * (wp * ws - wm * wq - 1) - 3 * wm * (n ** 2 * wm - (m1 - 1) ** 2 * wq)
        d[4] = (2 * m1 - 5) * (m1 - 2) * (wm - wq) + 3 * wm * (n ** 2 - (m1 - 2) ** 2)
        d[5] = n ** 2 - (m1 - 3) ** 2        
        tmp = 0
        for mu in range(0, 5):
            tmp = tmp + d[mu] * beta[m1 + 3 - (mu + 1)]
        beta[m1 - 3] = tmp / d[5]
    tmp = 0
    for m1 in range(0, n + 1):
        tmp = tmp + beta[m1]
    b = np.zeros((n + 1, 1))
    for m1 in range(0, n + 1):
        b[m1] = (-1) ** p * (beta[m1] / tmp)        
    return b

if __name__ == '__main__':        
    
    # parameters of the 'Zolotarev' polynomial
    n = 4; N = 2 * n + 1 # 'N' is the order of the polynomial
    R = 10 # 'R' is the value of the peak within [-1,+1] 
    
    # Getting the value of the parameter 'm' from a given 'R' (remember!, m=k^2)
    m = z_m_frm_R(N, R)
    print 'm =', m
    
    # Testing the side-lobe ratio (i.e., SLR in linear scale)
    R = z_Zolotarev_x2(N, m)
    print 'R =', R
    print 'SLR =', 10 * mp.log10(R ** 2), '(dB)' # SLR depends only on the magnitude of R here
    
    # x1, x2, x3 values ... just for the plotting purpose
    x1, x2, x3 = z_x123_frm_m(N, m)
    x11 = z_str2num(mp.nstr(x1, n=6))
    x22 = z_str2num(mp.nstr(x2, n=6))
    x33 = z_str2num(mp.nstr(x3, n=6))    
    
    # Evaluating the polynomial at some discrete points    
    x = np.linspace(-1.04, 1.04, num=500)
    y = []
    for i in range(len(x)):
        tmp = mp.nstr(z_Zolotarev(N, x[i], m), n=6)
        tmp = z_str2num(tmp)
        y.append(tmp)
    
    # Plotting the data obtained at those discrete points
    import matplotlib.pyplot as plt
    f1 = plt.figure(1)
    p1 = plt.subplot(111)
    p1.plot(x, y, '-r')
        
    # Getting the Zolotarev polynomial coefficients and roots
    coef, roots = z_Zolotarev_poly(N, m)
    print 'polynomial coefficients:', '\n', coef, '\n', 'polynomial roots:', '\n', roots    
    
    # Cross-checking the obtained coefficients
    x1 = np.linspace(-1.04, 1.04, num=100)
    y1 = np.polyval(coef, x1)
    
    # plotting the fitted polynomial
    p1.plot(x1, y1, '--k')
    p1.axhspan(-1, 1, facecolor='k', alpha=0.2)
    p1.axvspan(-x33, x33, facecolor='k', alpha=0.1)
    p1.axis('tight')
    p1.grid(True)
    plt.xlabel('$x$', fontsize=20)
    plt.ylabel(r'$Z$ $(x,$ $m)$', fontsize=20)
    p1.axis(fontsize=50)
    plt.show()    

#    # Testing the 'Zolotarev2' function 
#    p = 3; q = 6; n = p + q
#    k = mpf('0.682'); m = k ** 2
#    b = Zolotarev2(p, q, m)
#    print b
#    
#    x = np.arange(-1, 2, 0.01, dtype=float)
#    y = np.polyval(b, x)
#
#    import matplotlib.pyplot as plt            
#    f1 = plt.figure(1)
#    p1 = plt.subplot(111)    
#    p1.plot(x, y)
#    p1.axis('tight')
#    p1.grid(True)    
#    plt.xlabel(r'$x$')
#    plt.ylabel(r'$y$')
#    plt.show()