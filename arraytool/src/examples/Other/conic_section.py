#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

""" Program to draw a general conic section ...

for simple theoretical analysis refer:
http://zinka.files.wordpress.com/2010/06/conic-sections.pdf
"""

from __future__ import division
import matplotlib.pyplot as plt
from numpy import sqrt, arange, pi, cos, sin, set_printoptions, nan 

set_printoptions(precision=2, threshold=nan, edgeitems=None,
                 linewidth=None, suppress=None, nanstr=None, infstr=None)

#==============================================================================
# Parabola (Tested OK)
#==============================================================================

#a = 2;
#t = arange(-10, 10, 0.001)
#phi = arange(0, 2 * pi, 0.001)
#rho = (2 * a) / (1 - cos(phi))
#x1 = rho * cos(phi)
#y1 = rho * sin(phi)
#
#plt.figure(1)
#plt.plot(t, sqrt(4 * a * (t + a)), 'k',
#         t, -sqrt(4 * a * (t + a)), 'k',
#         x1, y1, '--r')         
#plt.axis([-20, 20, -20, 20])
#plt.grid(True)
#
#plt.show()

#==============================================================================
# Ellipse (Tested OK)
#==============================================================================

#a=6;
#b=4;
#t = arange(-10, 10, 0.001)
#phi = arange(0, 2*pi, 0.001)
#rho = (b**2)/(a-sqrt(a**2-b**2)*cos(phi))
#x1 = rho*cos(phi)
#y1 = rho*sin(phi)
#
#plt.figure(1)
##plt.plot(t, sqrt(b**2*(1-((t-sqrt(a**2-b**2))/(a))**2)),'k',
##         t, -sqrt(b**2*(1-((t-sqrt(a**2-b**2))/(a))**2)), 'k',
##         +sqrt(a**2-b**2), 0, 'o',
##         x1, y1, 'k')
#plt.plot(x1, y1, 'k',
#         +sqrt(a**2-b**2), 0, 'o',
#         +sqrt(a**2-b**2)-(a**2)/(sqrt(a**2-b**2)), 0, 'd')
#plt.axis([-20, 20, -20, 20])
#plt.grid(True)
#plt.show()

#==============================================================================
# Hyperbola (Tested OK)
#==============================================================================

a = 6;
b = 4;
t = arange(-34.5, 20, 0.001)
phi = arange(0, 2 * pi, 0.001)
rho = (b ** 2) / (a - sqrt(a ** 2 + b ** 2) * cos(phi))
x1 = rho * cos(phi)
y1 = rho * sin(phi)

plt.figure(1)
plt.plot(t, sqrt(b ** 2 * (-1 + ((t + sqrt(a ** 2 + b ** 2)) / (a)) ** 2)), 'k',
         t, -sqrt(b ** 2 * (-1 + ((t + sqrt(a ** 2 + b ** 2)) / (a)) ** 2)), 'k',
         - sqrt(a ** 2 + b ** 2), 0, 'o',
         - sqrt(a ** 2 + b ** 2) - (a ** 2) / (sqrt(a ** 2 + b ** 2)), 0, 'd',
         - sqrt(a ** 2 + b ** 2) + (a ** 2) / (sqrt(a ** 2 + b ** 2)), 0, 'd')
plt.plot(-sqrt(a ** 2 + b ** 2), 0, 'o',
         x1, y1, '--r')
plt.axis([-34.5, 20, -20, 20])
plt.grid(True)
plt.show()

