#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2014 Srinivasa Rao Zinka
# License: New BSD License.

import numpy as np
from scipy import integrate

def FS(fun_str_re, fun_str_im='0', T0=2 * np.pi, m_start= -5, m_stop=5, err_lim=1e-8):
    """Function to generate a finite number of Fourier series coefficients of 
    a periodic function."""
    
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

def IFS(FS, T0=2 * np.pi, m_start= -4, m_stop=4, x_min=0, x_max=2 * np.pi, x_num=10):
    """Function to reconstruct (or check) the periodic function from the 
    obtained Fourier coefficients"""
        
    m = np.arange(m_start, m_stop + 1)
    m = np.reshape(m, (-1, m.size))
    M = np.tile(m, (x_num, 1))
    x = np.linspace(x_min, x_max, num=x_num)
    x = np.reshape(x, (x.size, -1))
    X = np.tile(x, (1, m.size))
    FS = np.reshape(FS, (FS.size,-1))
    # Evaluating the inverse of the Fourier series
    IFS = np.dot(np.exp(1j * M * (2 * np.pi / T0) * X), FS)    
    
    return x, IFS

if __name__ == '__main__':
    
    