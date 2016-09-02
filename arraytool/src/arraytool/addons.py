#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2015 Srinivasa Rao Zinka
# License: New BSD License.

import planar as pl
import numpy as np
import matplotlib.pyplot as plt

def circ_bound (array_ip, radius = 10, plot=True):
    """ Function doc """

    for i in range(np.size(array_ip,0)):
        x = array_ip[i,0]
        y = array_ip[i,1]
        ele_dist = np.sqrt(x**2+y**2)
        if (ele_dist>radius):
            array_ip[i,3] = 0

    if  (plot):
        plt.figure(1)
        plt.scatter(array_ip[:,0], array_ip[:,1], c=array_ip[:,3], s=array_ip[:,3])
        plt.axis('equal')
        plt.grid(True)
        plt.show() # TODO draw circular boundary & show the count number in the plot itself

    ele_count = np.sum(array_ip[:,3])
    print ele_count

    return array_ip, ele_count

def quant ():
    """ Function doc """
    return

def beamwid (tp, F_dB, dB=-3):
    """ Function doc """

    F_peak = F_dB.max()
    msk1 = F_dB > F_peak+dB
    tp = tp*msk1
    minval = np.min(tp[np.nonzero(tp)])
    maxval = np.max(tp[np.nonzero(tp)])
    beamwid = maxval-minval

    return beamwid

def directivity_2d (tp, AF):
    """ Function doc """
    return

if __name__ == '__main__':

    # frequency and array-arrangement (actual values)
    freq = 15e9  # frequency of operation in Hzs
    wav_len = 3e8 / freq  # wavelength in meters
    M = 35  # no. of elements along the x-axis
    N = 35  # no. of elements along the y-axis
    a1 = 12e-3  # separation between elements along the x-axis in meters
    b1 = 12e-3  # separation between elements along the y-axis in meters
    gamma = np.pi / 3  # lattice angle in radians

    # normalized values
    a = a1 / wav_len  # 'a1' in-terms of lambda (wavelength)
    b = b1 / wav_len  # 'b1' in-terms of lambda (wavelength)

    #A = np.random.rand(N, M)
    A = np.ones((N, M))  # Uniform excitation
    array_ip = pl.ip_format(a, b, A, gamma, plot=False, stem=True, mayavi_app=False)

    circ_bound(array_ip, radius = 8.75, plot=True)

    # Calling the 'pattern_p' function to evaluate and plot 2D AF/GF/NF/NF0
    # tht, F = pl.pattern_p(array_ip, tht_scan=0, phi_scan=0, phi=0 * np.pi,
                   # tht_min=2, tht_max=4, tht_num=200, scale="dB",
                   # dB_limit=-40, factor="NF0", plot_type="rect", color='b',
                   # linewidth=1, linestyle='-', alpha=1, show=True)

    # print 180/np.pi*beamwid(tht, F)


# TODO amplitude & phase quantization (call the function as 'quant')
# (broadside array_ip is input; scan angles are specfied; A & P quantizations true or not;)
