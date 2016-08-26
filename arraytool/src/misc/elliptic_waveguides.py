#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2014 Srinivasa Rao Zinka
# License: New BSD License.

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import special as sp
from scipy import optimize

a = 111e-3
b = 74e-3

e = np.sqrt(1-b**2/a**2)
z = np.arccosh(1/e)

n = 0 # Order
#z = np.linspace(0, 2, num=100)
#z = 1.52
q = np.linspace(0, 7, num=50)

print e
print z
print q

Mc1,Mc11 = sp.mathieu_modcem1(n,q,z) # EVEN functions (modified or radial or CAPITAL)
Ms1,Ms11 = sp.mathieu_modsem1(n,q,z) # ODD functions (modified or radial or CAPITAL)

# Plotting TE modes
#plt.plot(q,np.sqrt(np.pi/2)*Mc11) # EVEN TE
#plt.plot(q,np.sqrt(np.pi/2)*Ms11) # ODD TE

# Plotting TM modes
#plt.plot(q,np.sqrt(np.pi/2)*Mc1) # EVEN TM
#plt.plot(q,np.sqrt(np.pi/2)*Ms1) # ODD TM

#Mcn1,Mcn11 = sp.mathieu_cem(n,q,z) # EVEN functions (normal or angular or SMALL)
#Msn1,Msn11 = sp.mathieu_sem(n,q,z) # ODD functions (normal or angular or SMALL)

# Plotting TE modes
#plt.plot(q,np.sqrt(np.pi/2)*Mcn11) # EVEN TE
#plt.plot(q,np.sqrt(np.pi/2)*Msn11) # ODD TE

plt.grid(True); plt.show()

# Obtaining roots ...
fun1 = lambda x: sp.mathieu_modcem1(n,x,z)[1] # EVEN TE
#fun1 = lambda x: sp.mathieu_modsem1(n,x,z)[1] # ODD TE
#fun1 = lambda x: sp.mathieu_modcem1(n,x,z)[0] # EVEN TM
#fun1 = lambda x: sp.mathieu_modsem1(n,x,z)[0] # ODD TM

print optimize.fsolve(fun1, 1.4)