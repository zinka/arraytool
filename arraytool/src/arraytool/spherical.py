#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2014 Srinivasa Rao Zinka
# License: New BSD License.

from __future__ import division
import numpy as np
from mayavi import mlab
from scipy import integrate
from scipy.special import sph_harm

# adjusting "matplotlib" label fonts
from matplotlib import rc
rc('text', usetex=True)

def SH(fun_str_re, fun_str_im='0', T0=2 * np.pi, m_start= -5, m_stop=5, err_lim=1e-8):
    r"""
    Function to generate a finite number of spherical harmonic series
    coefficients of a periodic function represented in (tht,phi) domain.
    """
    
    N = m_stop - m_start + 1
    FS = np.zeros((N, 1), dtype='complex')
    m_index = range(m_start, m_stop + 1)
    w0 = 2 * np.pi / T0

    for m in m_index:
        fun_re = lambda x: (eval(fun_str_re)) * np.cos(m * w0 * x) + (eval(fun_str_im)) * np.sin(m * w0 * x)
        fun_img = lambda x:-(eval(fun_str_re)) * np.sin(m * w0 * x) + (eval(fun_str_im)) * np.cos(m * w0 * x)       
        FS_re = integrate.quad(fun_re, 0, 2 * np.pi)
        FS_img = integrate.quad(fun_img, 0, 2 * np.pi)
        if ((FS_re[1] + FS_img[1]) < err_lim):
            FS[m - m_start] = (1 / T0) * (FS_re[0] + 1j * FS_img[0])
        else:
            print "Absolute error of the integration is not less than 1e-10 while calculating Fourier series"
            print "error(FS_re): ", FS_re[1]
            print "error(FS_img): ", FS_img[1]
        m_index = np.array(m_index) * (2 * np.pi / T0)
        m_index = np.reshape(m_index, (m_index.size, -1))
        
    return m_index, FS

if __name__ == '__main__':
    
    m = 2; n=5
    
    # Create a sphere
    r = 0.3
    pi = np.pi
    cos = np.cos
    sin = np.sin
    # tht and phi definitions are interchanged in scipy
    phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]
    
    # tht and phi definitions are interchanged in scipy
    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)
    
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
    mlab.clf()
    
    s = sph_harm(m, n, theta, phi).real
    mlab.mesh(x, y, z, scalars=s)
    mlab.axes(xlabel="X", ylabel="Y", zlabel="Z")
    mlab.show() 