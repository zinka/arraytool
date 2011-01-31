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