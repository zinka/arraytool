#! /usr/bin/env python

#def z_f(a, b, m):
#    """Summation of the series (eq.(23),[1])"""
#    K = ellipk(m)
#    K1 = ellipk(1 - m)
#    q1 = mp.exp(-mp.pi * K / K1)
#    NM = mp.cosh((a + b) * mp.pi / (2 * K))
#    DM = mp.cosh((a - b) * mp.pi / (2 * K))
#    fun = lambda r: ((-1) ** r / r) * ((q1 ** (2 * r)) / (1 - q1 ** (2 * r)))\
#     * mp.sinh(mp.pi * r * a / K1) ** mp.sinh(mp.pi * r * b / K1)
#    f = -(mp.pi * a * b / (K * K1)) + mp.log(NM / DM) - 4 * nsum(fun, [1, inf])
#    return f



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
