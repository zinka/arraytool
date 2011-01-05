#! /usr/bin/env python
""" Simple script for comparing different antenna polarization definitions """

from __future__ import division
import numpy as np
from pylab import *

np.set_printoptions(precision=2, edgeitems=None,
                 linewidth=None, suppress=None, nanstr=None, infstr=None)

#==============================================================================
# Basic parameters
#==============================================================================
phi_R = 0
e_act = 0
e_ref = 1
#==============================================================================
# Grid generation with "N" number of samples (denoted by the Nj)
#==============================================================================
[tht, phi] = np.mgrid[0:np.pi:100j, 0:2 * np.pi:100j]
#print tht, '\n', '\n', phi
#==============================================================================
# Actual "zeta" values
#==============================================================================
zeta_act = np.arctan (\
                      (np.tan(phi - phi_R) * (e_act + np.cos(tht)))\
                       / (1 + e_act * np.cos(tht))\
                       )
print zeta_act
#==============================================================================
# Reference "zeta" values
#==============================================================================
zeta_ref = np.arctan (\
                      (np.tan(phi - phi_R) * (e_ref + np.cos(tht)))\
                       / (1 + e_ref * np.cos(tht))\
                       )
#==============================================================================
# Directional error (zeta_act-zeta_ref)
#==============================================================================
dir_error = np.abs((zeta_act - zeta_ref))
#print dir_error
#==============================================================================
# Converting angles from "radians" to "degrees"
#==============================================================================
phi = phi* (180 / np.pi)
tht = tht* (180 / np.pi)
dir_error = dir_error* (180 / np.pi)
#==============================================================================
# Plotting the error as a contour plot
#==============================================================================
levels = np.linspace(0, 180, 10)
plt.figure(1)
plt.subplot(111)

CS = plt.contourf(phi, tht, dir_error, levels, cmap=cm.bone, alpha=0.5)
#CS = plt.contourf(phi, tht, dir_error, cmap=cm.bone, alpha=0.5)
cbar = colorbar(CS)
cbar.ax.set_ylabel('Directional Error')

plt.axis('tight')
#plt.title('Directional Error')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\theta$')

plt.show()