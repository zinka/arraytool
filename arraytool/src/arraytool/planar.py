#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

"""
Module for analysis and design of (uni-directional) planar phased array
antennas.

Progress:
Just some important basic routines are done. There is much more to be done!

References:
[] Will be provided later.
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from enthought.mayavi import mlab
import Zolotarev as Zol

def cutoff(F, dB_limit= -40):
    """
    when AF/GF/NF is 0, their dB value is '-infinity'. So, this function
    will be used to cut-off all the value below some 'dB_limit'.
    By default, 'dB_limit' is -40(dB).
    """
    msk1 = F < dB_limit
    fill = msk1 * dB_limit
    msk2 = F >= dB_limit
    F = F * (msk2) + fill
    return F

def ip_format(a, b, A, gamma=np.pi / 2, plot=False, color='b', linewidth=1,
              linestyle='-', alpha=1, show=True, stem=False, stemline='g--',
              stemmarker='ro', mayavi_app=False):
    """
    Function to generate the 'Arraytool' input format.

    a            : separation between elements along the x-axis in wavelengths
    b            : separation between elements along the y-axis in wavelengths
    A            : visual excitation matrix
    gamma        : lattice angle in radians
    plot         : if True, this function automatically produces a 2D/3D plot
                   showing the array excitation. In order to use this option
                   you need 'Matplotlib' and 'MayaVi' libraries.
    stem         : if True, a 'stem' representation of the excitation 
                   coefficients will be shown. By default, it is False.                   
    mayavi_app   : if True, the 3D plot will be opened in the MayaVi main
                   application itself.
    
    All other parameters are nothing but 'Matplotlib' parameters. These should
    be familiar to 'Matlab' or 'Matplotlib' users.                   
    """
    M = float(A.shape[1]) # no. of elements along the x-axis
    N = float(A.shape[0]) # no. of elements along the y-axis
    if (M == 1):  # i.e, linear array is along the y-direction
        a = 0
        gamma = np.pi / 2
    if (N == 1):  # i.e, linear array is along the x-direction
        b = 0
        gamma = np.pi / 2
    xlim = (M * a) / 2 # array is with in the x-limits [-xlim, +xlim]

    # Grid generation
    [x, y] = np.mgrid[0:M, 0:N]
    x = (x - (M - 1) / 2).T
    y = np.flipud((y - (N - 1) / 2).T)
    x = x * a
    y = y * b # rectangular grid is generated

    # modifying the rect-grid according to the given lattice angle 'gamma'
    if (gamma != np.pi / 2):
        x = x + (y / np.tan(gamma)) % a

    # Adjusting the rows so that the array lies within [-xlim, +xlim]
    for i1 in range(int(N)):
        if (x[i1, 0] < -xlim):
            x[i1, :] = x[i1, :] + a # RIGHT shifting the row by 'a'
        if (x[i1, -1] > xlim):
            x[i1, :] = x[i1, :] - a # LEFT shifting the row by 'a'

    # Finally, arranging all the data into 'Arraytool' input format
    x = np.reshape(x, (M * N, -1))
    y = np.reshape(y, (M * N, -1))
    z = np.zeros_like(x) # because only planar arrays are permitted here
    A = np.reshape(A, (M * N, -1))
    array_ip = np.hstack((x, y, z, A))  # finally, 'Arraytool' input format

    # plotting the 'absolute' value of the array excitation (2D/3D)
    if (plot):
        # checking whether 'A' has any imaginary values
        if((A.imag > 1e-10).sum()):
            A_plt = abs(A) # if A.imag are significant, then '|A|' will be plotted
        else:
            A_plt = A.real # if A.imag are negligible, then 'A'  will be plotted
        if (M == 1):  # i.e, linear array is along the y-direction
            plt.plot(y, A_plt, color=color, linewidth=linewidth,
                         linestyle=linestyle, alpha=alpha)
            if(stem): plt.stem(y, A_plt, linefmt=stemline, markerfmt=stemmarker)
            plt.axis('tight'); plt.grid(True)
            plt.xlabel(r'$y$'); plt.ylabel(r'$\left|A_{n}\right|$')
            if(show): plt.show()
        elif (N == 1):  # i.e, linear array is along the x-direction
            plt.plot(x, A_plt, color=color, linewidth=linewidth,
                         linestyle=linestyle, alpha=alpha)
            if(stem): plt.stem(x, A_plt, linefmt=stemline, markerfmt=stemmarker)
            plt.axis('tight'); plt.grid(True)
            plt.xlabel(r'$x$'); plt.ylabel(r'$\left|A_{m}\right|$')
            if(show): plt.show()
        else:
            if (mayavi_app): # this option opens the 3D plot in MayaVi Application
                mlab.options.backend = 'envisage'
            s1 = mlab.quiver3d(x, y, z, z, z, A_plt) # stem3D representation
            ranges1 = [x.min(), x.max(), y.min(), y.max(), A_plt.min(), A_plt.max()]
            mlab.axes(xlabel="x", ylabel="y", zlabel="Excitation", ranges=ranges1)
            s1.scene.isometric_view()
            if(show): mlab.show()
    return array_ip

def ATE(array_ip):
    """A simple function to evaluate the array taper efficiency (ATE)."""
    A = abs(array_ip[:, 3])
    A2 = A * A
    ATE = (A.sum())**2 / (len(A) * A2.sum())
    return ATE

def K_norm(array_ip):
    """A simple function to evaluate the normalized bore-sight slope."""
    X = array_ip[:, 0]
    A = array_ip[:, 3]
    AX = A * X
    A1 = abs(A)
    A2 = A1 * A1    
    K_norm = abs(AX.sum()) / np.sqrt(A2.sum())
    return K_norm

def AF_zeros(a, M, R, dist_type, nbar=False, alpha=0):
    """
    This function gives array-factor zeros corresponding to different
    types of array distributions.

    a         : separation between the elements along the x-axis in wavelengths
    M         : number of elements along the x-axis
    R         : side-lobe ratio in linear scale
    type      : type of the distribution, e.g., 'Dolph' for Dolph-Chebyshev
    nbar      : dilation parameter
    alpha     : Taylor's asymptotic tapering parameter (this 'alpha' has nothing
                to do with 'plotting alpha', i.e., the transparency parameter)
    """
    k = 2 * np.pi # (angular) wave-number, which is 2*pi when lambda = 1
    m = np.ceil((M - 2) / 2)
    n = np.arange(1, 1 + m, 1) # number of zeros for symmetric array-factors
    na = np.arange(1, M, 1) # number of zeros for 'asymmetric' array-factors
    
    if(dist_type == "Dolph"): # Dolph zeros
        c = np.cosh(np.arccosh(R) / (M - 1))
        U0 = (2 / (a * k)) * np.arccos((np.cos(np.pi * (2 * n - 1) / (2 * M - 2))) / c)        
    elif(dist_type == "Riblet"): # Riblet zeros
        c1 = np.cosh(np.arccosh(R) / m)
        c = np.sqrt((1 + c1) / (2 + (c1 - 1) * np.cos(k * a / 2) ** 2))
        alph = c * np.cos(k * a / 2)
        xi = (1 / c) * np.sqrt(((1 + alph ** 2) / 2) + ((1 - alph ** 2) / 2) * 
                               np.cos(((2 * n - 1) * np.pi) / (2 * m)))
        U0 = (2 / (a * k)) * np.arccos(xi)
    elif(dist_type == "DuhamelB"): # Duhamel bi-directional end-fire array zeros
        if(a < 0.5): c = np.cosh(np.arccosh(R) / (M - 1)) / np.sin((k * a) / 2)
        else: c = np.cosh(np.arccosh(R) / (M - 1))
        U0 = (2 / (a * k)) * np.arcsin((np.cos(np.pi * (2 * n - 1) / (2 * M - 2))) / c)
    elif(dist_type == "DuhamelU"): # Duhamel uni-directional end-fire array zeros
        Lamb = np.cosh(np.arccosh(R) / (M - 1))
        xi = (2 / a) * (0.5 * np.pi - (np.arctan(np.tan(k * a / 2) * ((Lamb + 1) / (Lamb - 1)))))
        c = 1 / (np.sin((xi - k) * a / 2))
        U0 = -(xi / k) + (2 / (a * k)) * np.arcsin((
                            np.cos(np.pi * (2 * na - 1) / (2 * M - 2))) / c)
    elif(dist_type == "McNamara-s"): # McNamara-Zolotarev sum-pattern zeros
        U0 = "Yet to be done"
    elif(dist_type == "McNamara-d"): # McNamara-Zolotarev difference-pattern zeros
        if(a < 0.5): c = 1 / np.sin((k * a) / 2)
        else: c = 1
        m1 = Zol.z_m_frm_R(M - 1, R)
        xn = Zol.z_Zolotarev_poly(N=M - 1, m=m1)[1][m + 1:]
        U0 = (2 / (a * k)) * np.arcsin(xn / c)
    if(nbar): # Taylor's Dilation procedure
        # see if you can change the below LONG logic
        if((dist_type == "Dolph") or (dist_type == "Riblet") or (dist_type == "McNamara-s")):
            n_gen = np.arange(nbar, 1 + m, 1) # indices of the generic zeros
            U0_gen = (n_gen + alpha / 2) * (1 / (M * a)) # generic sum zeros
        elif(dist_type == "McNamara-d"):
            # THIS NEEDS MODIFICATION ... some thing wrong!
            n_gen = np.arange(nbar, 1 + m, 1) # indices of the generic zeros
            U0_gen = (n_gen + (alpha + 1) / 2) * (1 / (M * a)) # generic difference zeros
        sigma = U0_gen[0] / U0[nbar - 1] # Dilation factor
        U0 = np.hstack((sigma * U0[0:nbar - 1], U0_gen)) # Dilated zeros
    U0 = np.reshape(U0, (len(U0), -1))
    return U0

def A_frm_zeros(U0, a, M, symmetry):
    """
    This function gives array excitation coefficients corresponding to the
    given array factor zeros.

    U0        : arrayfactor zeros... usually in the format of the output of the
                'AF_zeros' function.
    a         : separation between elements along the x-axis in wavelengths
    M         : number of elements along the x-axis
    symmetry  : whether the array excitation is symmetric or not. this simplifies
                the numerical evaluation process. Allowed values are 'even',
                'odd' and False.
    """
    k = 2 * np.pi # (angular) wave-number, which is 2*pi when lambda = 1
    sz = len(U0)
    UU = np.tile(U0, sz)
    if(symmetry == "even"):
        if(M % 2 == 0):
            tmp1 = np.arange(1, 1 + sz, 1) + 0.5
            tmp2 = 2 * np.cos(k * U0 * a / 2)
        else:
            tmp1 = np.arange(1, 1 + sz, 1)
            tmp2 = np.ones_like(U0)
        tmp1 = np.reshape(tmp1, (-1, sz))
        TT = np.tile(tmp1, (sz, 1))
        CC = -np.linalg.inv(2 * np.cos(k * UU * TT * a))
        A = np.dot(CC, tmp2)
        A1 = np.flipud(A)
        if(M % 2 == 0):
            A_tot = np.vstack((A1, 1, 1, A))
        else:
            A_tot = np.vstack((A1, 1, A))
    elif(symmetry == "odd"):
        if(M % 2 == 0):
            tmp1 = np.arange(1, 1 + sz, 1) + 0.5
            tmp2 = 2 * np.sin(k * U0 * a / 2)
        else:
            tmp1 = np.arange(1, 1 + sz, 1)
            tmp2 = np.ones_like(U0)
        tmp1 = np.reshape(tmp1, (-1, sz))
        TT = np.tile(tmp1, (sz, 1))
        CC = -np.linalg.inv(2 * np.sin(k * UU * TT * a))
        A = np.dot(CC, tmp2)
        A1 = -np.flipud(A)
        if(M % 2 == 0):
            A_tot = np.vstack((A1, -1, 1, A))
        else:
            A_tot = np.vstack((A1, 1, A))
    elif(symmetry == False):
        tmp1 = np.arange(1, 1 + sz, 1) - (M - 1) / 2
        tmp1 = np.reshape(tmp1, (-1, sz))
        TT = np.tile(tmp1, (sz, 1))
        CC = -np.linalg.inv(np.exp(1j * k * UU * TT * a))
        tmp2 = np.exp(-1j * k * U0 * ((M - 1) / 2) * a)
        A = np.dot(CC, tmp2)
        A_tot = np.vstack((1, A))
    return A_tot

def pattern_u(array_ip, u_scan=0, u_min= -1, u_max=1, u_num=50, scale="dB",
              dB_limit= -40, factor="GF", plot_type="rect", lattice=False,
              color='b', linewidth=1, linestyle='-', alpha=1, show=True):
    """
    Function to evaluate 2d array-factor(AF) or gain-factor(GF) of a
    linear array in u-domain. By default, this function calculates gain-factor.

    array_ip      : input data in 'Arraytool' input format
    u_scan        : beam scan position. For example, if you need to scan to
                    theta=30deg, u_scan = sin(pi/6). By default, 'u_scan' is 0.
    u_min, u_max  : limits of u-domain, by default they are -1 and +1 which
                    correspond to the visible-space
    u_num         : number of sampling points between 'u_min' and 'u_max',
                    including boundaries. By default, is 50.
    scale         : By default, this is "dB". If 'scale="linear"', this function
                    produces a 2D plot in linear scale.
    dB_limit      : when AF/GF/NF is 0, their dB value is -infinity. So, we will
                    cut-off all the value below some 'dB_limit'. By default,
                    'dB_limit' is -40(dB).
    factor        : default value is gain-factor, "GF". If you want
                    array-factor or normalized factor, choose "AF" or "NF".
    plot_type     : by default, is "rect". If you want polar plot, use
                    'plot_type="polar"'. Finally, if 'plot_type' is False, it
                    doesn't plot anything. You need 'matplotlib' library to use
                    this option.
    lattice       : by default, is False. If True, it will highlight both
                    visible-space and lattice period in u-domain of
                    "rectangular pattern" plot mode. Lattice period has meaning
                    only if the array is "uniformly space".
                    
    All other parameters are nothing but 'Matplotlib' parameters. These should
    be familiar to 'Matlab' or 'Matplotlib' users.
    """
    x = array_ip[:, 0]
    y = array_ip[:, 1]
    z = array_ip[:, 2]
    A = array_ip[:, 3] # un-packing "array_ip" finished
    k = 2 * np.pi # (angular) wave-number, which is 2*pi when lambda = 1

    # Making sure all elements in y and z columns of the "array_ip" are zeros
    z_flag = True
    if ((abs(y) > 0).sum()):
        print "All elements in y-column of array input should be zero."
        z_flag = False
    elif ((abs(z) > 0).sum()):
        print "All elements in z-column of array input should be zero."
        z_flag = False

    # After making sure, proceed to the next level, i.e., evaluate the pattern
    if(z_flag):

        u = np.linspace(u_min, u_max, num=u_num)
        u = np.reshape(u, (u_num, -1))
        A = np.reshape(A, (len(A), -1))
        U = np.tile(u - u_scan, len(x))
        X = np.tile(x, (u_num, 1))

        # Evaluating array-factor of the linear array -2
        AF = np.dot(np.exp(1j * k * U * X), A)

        # Evaluation of F = (AF/GF/NF) => depending upon the user's choice
        if(factor == "AF"):
            F = AF; n1 = ""; ff = "Array-Factor "; f1 = "AF "
        elif(factor == "GF"):
            P_inc = ((abs(A)) ** 2).sum()
            GF = AF / np.sqrt(P_inc) # Converting the AF to GF
            F = GF; n1 = ""; ff = "Gain-Factor "; f1 = "GF "
        elif(factor == "NF"):
            norm_fact = (abs(A)).sum()
            F = AF / norm_fact
            n1 = "Normalized "; ff = "Factor "; f1 = "NF "

        # converting 'F' from linear to dB scale, if needed
        if(scale == "linear"):
            F_plt = abs(F)
            ss = "in linear scale"
        elif(scale == "dB"):
            F = 20 * np.log10(abs(F))
            # cutoff the "F" below some limit ... just for the plotting purpose
            F_plt = cutoff(F, dB_limit)
            ss = "in dB scale"

        # plotting the factor 'F_plt'
        if(plot_type):
            if(plot_type == "rect"): # rectangular plot
                plt.plot(u, F_plt, color=color, linewidth=linewidth,
                         linestyle=linestyle, alpha=alpha)
                if(lattice): # highlighting the visible-space and unit-lattice
                    plt.axvspan(-1, +1, facecolor='y', alpha=0.2)
                    lim = -np.pi / ((x[2] - x[1]) * k)
                    plt.axvspan(-lim, +lim, facecolor='b', alpha=0.2)
                plt.axis('tight'); plt.grid(True)
                plt.xlabel('u, where "u=sin(theta)" in visible-space')
                plt.ylabel(f1 + '(u)')
            if(plot_type == "polar"): # polar plot
                th = np.arcsin(u)
                if(scale == "linear"):
                    plt.polar(th, F_plt, color=color, linewidth=linewidth,
                              linestyle=linestyle, alpha=alpha)
                    plt.polar(np.pi - th, F_plt, color=color, linewidth=linewidth,
                              linestyle=linestyle, alpha=alpha)
                if(scale == "dB"):
                    plt.polar(th, F_plt - dB_limit, color=color,
                              linewidth=linewidth, linestyle=linestyle, alpha=alpha)
                    plt.polar(np.pi - th, F_plt - dB_limit, color=color,
                              linewidth=linewidth, linestyle=linestyle, alpha=alpha)
            plt.title(n1 + ff + ss)
            if(show): plt.show()
    return u, F

def pattern_uv(array_ip, u_scan=0, v_scan=0, u_min= -1, u_max=1, u_num=50,
               v_min= -1, v_max=1, v_num=50, uv_abs=1, scale="dB",
               dB_limit= -40, factor="GF", plot_type="rect", lattice=False,
               mayavi_app=False):
    """ Not finished yet. """
    x = array_ip[:, 0]
    y = array_ip[:, 1]
    z = array_ip[:, 2]
    A = array_ip[:, 3] # un-packing "array_ip" finished
    k = 2 * np.pi # (angular) wave-number, which is 2*pi when lambda = 1
    u_numj = complex(0, u_num)
    v_numj = complex(0, v_num)

    # Making sure all elements in the z-column of the "array_ip" are zeros
    z_flag = True
    if ((abs(z) > 0).sum()):
        print "All elements in the z-column of array input should be zero."
        z_flag = False

    # After making sure, proceed to the next level, i.e., evaluate the pattern
    if(z_flag):

        [u, v] = np.mgrid[u_min:u_max:u_numj, v_min:v_max:v_numj]
        u1 = np.reshape(u, (u.size, -1))
        v1 = np.reshape(v, (v.size, -1))
        A = np.reshape(A, (len(A), -1))
        U = np.tile(u1 - u_scan, len(x))
        V = np.tile(v1 - v_scan, len(x))
        X = np.tile(x, (u.size, 1))
        Y = np.tile(y, (u.size, 1))

        # Evaluating array-factor of the planar array
        AF1 = np.dot(np.exp(1j * k * (U * X + V * Y)), A)
        AF = np.reshape(AF1, u.shape)

        # Evaluation of F = (AF/GF/NF) => depending upon the user's choice
        if(factor == "AF"):
            F = AF; n1 = ""; ff = "Array-Factor "; f1 = "AF "
        elif(factor == "GF"):
            P_inc = ((abs(A)) ** 2).sum()
            GF = AF / np.sqrt(P_inc) # Converting the AF to GF
            F = GF; n1 = ""; ff = "Gain-Factor "; f1 = "GF "
        elif(factor == "NF"):
            norm_fact = (abs(A)).sum()
            F = AF / norm_fact
            n1 = "Normalized "; ff = "Factor "; f1 = "NF "

        # converting 'F' from linear to dB scale, if needed
        if(scale == "linear"):
            F_plt = abs(F)
            ss = "in linear scale"
        elif(scale == "dB"):
            F = 20 * np.log10(abs(F))
            # cutoff the "F" below some limit
            F = cutoff(F, dB_limit)
            F_plt = F
            ss = "in dB scale"

        # plotting the factor (AF/GF/NF)
        if(plot_type):
            if (mayavi_app): # opens the 3D plot in MayaVi Application
                mlab.options.backend = 'envisage'
            if(plot_type == "rect"): # rectangular plot
                me1 = mlab.surf(u, v, F_plt, warp_scale='auto')
                ranges1 = [u_min, u_max, v_min, v_max, F_plt.min(), F_plt.max()]
                mlab.axes(xlabel='u', ylabel='v', zlabel=f1 + '(u,v)',
                          ranges=ranges1)
                me1.scene.isometric_view()
                mlab.show()
            if(plot_type == "polar"): # polar plot
                print "to be done."
    return x

if __name__ == '__main__':

    # frequency and array-arrangement (actual values)
    freq = 10e9 # frequency of operation in Hzs
    wav_len = 3e8 / freq # wavelength in meters
    M = 10 # no. of elements along the x-axis
#    N = 5 # no. of elements along the y-axis
    a1 = 17e-3 # separation between elements along the x-axis in meters
    b1 = 17e-3 # separation between elements along the y-axis in meters
    gamma = np.pi / 2 # lattice angle in radians

    # normalized values
    a = a1 / wav_len # 'a1' in-terms of lambda (wavelength)
    b = b1 / wav_len # 'b1' in-terms of lambda (wavelength)

    SLR = 25 # side-lobe ratio in dB
    R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale

#    # Array excitation coefficient matrix
#
#    A = np.array([1,2,3,4,5]) # manually entering the coefficients
#    A = np.reshape(A, (len(A),-1)).T
#
#    A = np.ones((N, M)) # Uniform excitation
#
#    A = np.random.rand(N, M) # Random excitation

    # Using the function 'AF_zeros' to find arrayfactor zeros
    U0 = AF_zeros(a, M, R, dist_type="McNamara-d", nbar=False, alpha=0)
    print 'arrayfactor zeros:', '\n', U0

    # Obtaining array excitation coefficients from the arrayfactor zeros
    A = A_frm_zeros(U0, a, M, symmetry="odd").T # finding excitation coefficients
    print 'array coefficients:', '\n', A.T

    # Converting the 'excitation & position' info into 'Arraytool' input format
    array_ip = ip_format(a, b, A, gamma, plot=False, stem=True, mayavi_app=False)
    ATE = ATE(array_ip) # array taper efficiency

#    # Calling the 'pattern_u' function to evaluate and plot 2D AF/GF/NF
#    pattern_u(array_ip, u_scan=0, u_min= -1, u_max=1, u_num=500, scale="dB",
#              dB_limit= -60, factor="AF", plot_type="rect", lattice=True)

    # Calling the 'pattern_uv' function to evaluate and plot 3D AF/GF/NF
    pattern_uv(array_ip, u_scan=0, v_scan=0, u_min= -1, u_max=1, u_num=100,
               v_min= -1, v_max=1, v_num=100, uv_abs=1, scale="dB",
               dB_limit= -40, factor="NF", plot_type="rect", lattice=False,
               mayavi_app=False)

#==============================================================================
# Programming tasks
#==============================================================================
# odd symmetry ... optimize A_from_zeros
# use backslah instead of inverse

