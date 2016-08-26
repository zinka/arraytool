#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2014 Srinivasa Rao Zinka
# License: New BSD License.

import numpy as np
import planar as pl
import matplotlib.pyplot as plt

M = 9  # no. of elements along the x-axis
N = 1  # no. of elements along the y-axis
gamma = np.pi / 2  # lattice angle in radians

# normalized separation values
a = 0.5  # 'a1' in-terms of lambda (wavelength)
b = 0.5  # 'b1' in-terms of lambda (wavelength)

A = np.ones((N, M))  # Uniform excitation  
array_ip = pl.ip_format(a, b, A, gamma, plot=False, stem=True, mayavi_app=False)

#==============================================================================
# Orthogonal beams
#==============================================================================

u1, F1 = pl.pattern_u(array_ip, u_scan=0, u_min=-3, u_max=3, u_num=700, scale="linear",
          dB_limit=-40, factor="NF0", plot_type=False, lattice=True)

u2, F2 = pl.pattern_u(array_ip, u_scan=1/(M*a), u_min=-3, u_max=3, u_num=700, scale="linear",
          dB_limit=-40, factor="NF0", plot_type=False, lattice=True)

u3, F3 = pl.pattern_u(array_ip, u_scan=2/(M*a), u_min=-1, u_max=1, u_num=700, scale="linear",
          dB_limit=-40, factor="NF0", plot_type=False, lattice=True)

u4, F4 = pl.pattern_u(array_ip, u_scan=3/(M*a), u_min=-1, u_max=1, u_num=700, scale="linear",
          dB_limit=-40, factor="NF0", plot_type=False, lattice=True)

u5, F5 = pl.pattern_u(array_ip, u_scan=4/(M*a), u_min=-1, u_max=1, u_num=700, scale="linear",
          dB_limit=-40, factor="NF0", plot_type=False, lattice=True)

u6, F6 = pl.pattern_u(array_ip, u_scan=5/(M*a), u_min=-1, u_max=1, u_num=700, scale="linear",
          dB_limit=-40, factor="NF0", plot_type=False, lattice=True)

u7, F7 = pl.pattern_u(array_ip, u_scan=6/(M*a), u_min=-1, u_max=1, u_num=700, scale="linear",
          dB_limit=-40, factor="NF0", plot_type=False, lattice=True)

u8, F8 = pl.pattern_u(array_ip, u_scan=7/(M*a), u_min=-1, u_max=1, u_num=700, scale="linear",
          dB_limit=-40, factor="NF0", plot_type=False, lattice=True)

u9, F9 = pl.pattern_u(array_ip, u_scan=8/(M*a), u_min=-3, u_max=3, u_num=700, scale="linear",
          dB_limit=-40, factor="NF0", plot_type=False, lattice=True)

plt.plot(u1, F1)
plt.plot(u2, F2)
# plt.plot(u3, F3)
# plt.plot(u4, F4)
# plt.plot(u5, F5)
# plt.plot(u6, F6)
# plt.plot(u7, F7)
# plt.plot(u8, F8)
plt.plot(u9, F9)

plt.axis('tight')
plt.grid(True)
plt.xlabel(r'$\sin \theta$')
plt.ylabel('AF (Linear Scale)')
plt.show()

#==============================================================================
# Beam broadening
#==============================================================================

# u1, F1 = pl.pattern_u(array_ip, u_scan=0, u_min=-1, u_max=1, u_num=700, scale="linear",
#           dB_limit=-40, factor="NF0", plot_type=False, lattice=True)
# 
# # u2, F2 = pl.pattern_u(array_ip, u_scan=0.5, u_min=-1, u_max=1, u_num=700, scale="linear",
# #           dB_limit=-40, factor="NF0", plot_type=False, lattice=True)
# 
# u3, F3 = pl.pattern_u(array_ip, u_scan=1/np.sqrt(2), u_min=-1, u_max=1, u_num=700, scale="linear",
#           dB_limit=-40, factor="NF0", plot_type=False, lattice=True)
# 
# plt.plot(u1, abs(F1))
# # plt.plot(u2, abs(F2))
# plt.plot(u3, abs(F3))
# 
# plt.axis('tight')
# plt.grid(True)
# plt.xlabel(r'$\sin \theta$')
# plt.ylabel('AF (Linear Scale)')
# plt.show()
# 
# pl.pattern_p(array_ip, tht_scan=(np.pi/6), phi_scan=0, phi=0 * np.pi,
#                    tht_min=0, tht_max=2 * np.pi, tht_num=200, scale="dB",
#                    dB_limit=-40, factor="NF0", plot_type="polar", color='b',
#                    linewidth=1, linestyle='-', alpha=1, show=True)