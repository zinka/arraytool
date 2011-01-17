#! /usr/bin/env python
""" Simple script for comparing different antenna polarization definitions """

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from enthought.mayavi import mlab

np.set_printoptions(precision=2, edgeitems=None,
                 linewidth=None, suppress=None, nanstr=None, infstr=None)

#==============================================================================
# Basic parameters
#==============================================================================
e1 = 1
e2 = 1000000
phi_R_1 = 0
phi_R_2 = 0
#==============================================================================
# Grid generation with "N" number of samples (denoted by the Nj)
#==============================================================================
[tht, phi] = np.mgrid[0:(1 * np.pi):100j, 0:(2 * np.pi):100j]
#==============================================================================
# sin(zeta1), cos(zeta1), sin(zeta2) and cos(zeta2) values
#==============================================================================
sc1 = (np.sin(phi - phi_R_1) * (np.cos(tht) + e1))
cc1 = (np.cos(phi - phi_R_1) * (1 + e1 * np.cos(tht)))
s1 = sc1 / np.sqrt(sc1 ** 2 + cc1 ** 2)
c1 = cc1 / np.sqrt(sc1 ** 2 + cc1 ** 2)

sc2 = (np.sin(phi - phi_R_2) * (np.cos(tht) + e2))
cc2 = (np.cos(phi - phi_R_2) * (1 + e2 * np.cos(tht)))
s2 = sc2 / np.sqrt(sc2 ** 2 + cc2 ** 2)
c2 = cc2 / np.sqrt(sc2 ** 2 + cc2 ** 2)
#==============================================================================
# sin(zeta1-zeta2) and cos(zeta1-zeta2) values
#==============================================================================
s12 = s1 * c2 - s2 * c1
c12 = c1 * c2 + s1 * s2
#==============================================================================
# (zeta1-zeta2) value: 0 <= (zeta1-zeta2) <= 2pi
#==============================================================================
z12 = ((s12 > 0) * np.arccos(c12)) + ((s12 < 0) * (2 * np.pi - np.arccos(c12)))
#==============================================================================
# (zeta1-zeta2) value: 0 <= (zeta1-zeta2) <= pi
#==============================================================================
#z12 = np.arccos(c12)
#==============================================================================
# (zeta1-zeta2) value: 0 <= (zeta1-zeta2) <= pi/2
#==============================================================================
#z12 = np.arccos(c12) # modify this code later
#==============================================================================
# making the figure smooth by removing the "Nan" values
#==============================================================================
#==============================================================================
# Converting angles from "radians" to "degrees"
#==============================================================================
phi = phi * (180/np.pi)
tht = tht * (180/np.pi)
z12 = z12 * (180/np.pi)
#==============================================================================
# Plotting "(zeta1-zeta2)" as a contour plot
#==============================================================================
##fi1 = plt.figure(1)
#pl1 = plt.subplot(111)
#
#levels = np.linspace(0, 360, 19)
#
## Contour with filled area
#cpf1 = pl1.contourf(tht, phi, z12, levels, cmap = plt.cm.bone,  alpha=0.5)
#cbar = plt.colorbar(cpf1)
#cbar.ax.set_ylabel(r'$\Delta\zeta=\zeta_1-\zeta_2$ (degrees)')
#
### Contour without filled area
##plt.rcParams['contour.negative_linestyle'] = 'solid'
##cp1 = plt.contour(tht, phi, z12, levels, colors='k', alpha=1, linewidths=1)
##plt.clabel(cp1, fontsize=9, inline=True, fmt='%1.0f', manual=False)
#
##pl1.axis([0, 360, 0, 180])
#pl1.axis('tight')
#pl1.grid(True)
#
##plt.title('Directional Error (degrees)')
#plt.xlabel(r'$\theta$ (degrees)')
#plt.ylabel(r'$\phi$ (degrees)')
#plt.show()
#==============================================================================
# Plotting "(zeta1-zeta2)" as a 3D surface plot
#==============================================================================
fi1 = mlab.figure()
me1 = mlab.mesh(tht, phi, z12)
a1 = mlab.axes(xlabel="Theta", ylabel="Phi", zlabel="Error")
mlab.show()