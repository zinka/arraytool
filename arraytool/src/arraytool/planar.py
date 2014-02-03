#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2013 Srinivasa Rao Zinka
# License: New BSD License.

"""
Module for analysis and design of (uni-directional) planar phased array
antennas.

**Progress:**
Just some important basic routines are done. There is much more to be done!

**References**

- Will be provided later
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import Zolotarev as Zol
import Tkinter as ti
import tkFileDialog as tkdlg

# adjusting "matplotlib" label fonts ... can't save .svz files using this option
from matplotlib import rc
rc('text', usetex=True)

def at_import(dtype='complex'):
    r"""
    A simple function to import a CSV text file as a Numpy ndarray
    
    :param dtype: Data-type of the resulting array; default:complex. For further
                  information, see numpy.loadtxt.
                  
    :rtype:       ip, a Numpy ndarray
    """
    master = ti.Tk(); master.withdraw()  # hiding tkinter window 
    file_path = tkdlg.askopenfilename(title="Open file", filetypes=[("txt file",
                ".csv"), ("All files", ".*")]); master.quit()
    ip = np.loadtxt(file_path, delimiter=',', dtype=dtype)
    return ip

def at_export(data=np.ones((0, 0)), data_ID=False, fmt='%.4e', mode='a'):
    r"""
    A simple function to export a Numpy ndarray as a CSV text file. You can use
    this function to export "any Numpy array data" for further manipulations.
    
    :param data:        Numpy ndarray
    :param data_ID:     a string to represent the data being exported
    :param fmt:         str or sequence of strs. For more information, see
                        numpy.savetxt.
    :param mode:        file opening mode, e.g., 'w', 'a', etc
    
    :rtype:             A CSV text file 
    """
    master = ti.Tk(); master.withdraw()  # hiding tkinter window 
    file_path = tkdlg.asksaveasfile(mode, title="Save file", filetypes=[("txt file",
                ".csv"), ("All files", ".*")]); master.quit()
    if(data_ID):
        file_path.write(data_ID + '\n' + '\n')
    np.savetxt(file_path, data, delimiter=',', fmt=fmt)
    file_path.write('\n')
    return

def ip_format(a, b, A, gamma=np.pi / 2, plot=False, color='b', linewidth=1,
              linestyle='-', alpha=1, show=True, stem=False, stemline='g--',
              stemmarker='ro', fgcolor=(1, 1, 1), bgcolor=(0.5, 0.5, 0.5),
              mayavi_app=False):
    r"""
    Function to generate the standard 'Arraytool' input format.

    :param a:          separation between elements along the x-axis in wavelengths
    :param b:          separation between elements along the y-axis in wavelengths
    :param A:          visual excitation matrix
    :param gamma:      lattice angle in radians
    :param plot:       if True, produces a 2D/3D plot of the array excitation
    :param stem:       if True, the array excitation is plotted as 'stem plot'
    :param mayavi_app: if True, the 3D plot will be opened in the MayaVi application
                       so that we can save the visualization for future reference
    
    All other parameters are nothing but the 'Matplotlib' parameters. These
    should be familiar to 'Matlab' or 'Matplotlib' users.
    
    :rtype:            array_ip, a Numpy array of size (Number_elements(A),4)
    """
    M = float(A.shape[1])  # no. of elements along the x-axis
    N = float(A.shape[0])  # no. of elements along the y-axis
    if (M == 1):  # i.e, linear array is along the y-direction
        a = 0
        gamma = np.pi / 2
    if (N == 1):  # i.e, linear array is along the x-direction
        b = 0
        gamma = np.pi / 2
    xlim = (M * a) / 2  # array is with in the x-limits [-xlim, +xlim]

    # Grid generation
    [x, y] = np.mgrid[0:M, 0:N]
    x = (x - (M - 1) / 2).T
    y = np.flipud((y - (N - 1) / 2).T)
    x = x * a
    y = y * b  # rectangular grid is generated

    # modifying the rect-grid according to the given lattice angle 'gamma'
    if (gamma != np.pi / 2):
        x = x + (y / np.tan(gamma)) % a

    # Adjusting the rows so that the array lies within [-xlim, +xlim]
    for i1 in range(int(N)):
        if (x[i1, 0] < -xlim):
            x[i1, :] = x[i1, :] + a  # RIGHT shifting the row by 'a'
        if (x[i1, -1] > xlim):
            x[i1, :] = x[i1, :] - a  # LEFT shifting the row by 'a'

    # Finally, arranging all the data into 'Arraytool' input format
    x = np.reshape(x, (M * N, -1))
    y = np.reshape(y, (M * N, -1))
    z = np.zeros_like(x)  # because only planar arrays are permitted here
    A = np.reshape(A, (M * N, -1))
    array_ip = np.hstack((x, y, z, A))  # finally, 'Arraytool' input format

    # plotting the 'absolute' value of the array excitation (2D/3D)
    if (plot):
        # checking whether 'A' has any imaginary values
        if((abs(A.imag) > 1e-10).sum()):
            A_plt = abs(A)  # if A.imag are significant, then '|A|' will be plotted
        else:
            A_plt = A.real  # if A.imag are negligible, then 'A'  will be plotted
        if (M == 1):  # i.e, linear array is along the y-direction
            plt.plot(y, A_plt, color=color, linewidth=linewidth,
                         linestyle=linestyle, alpha=alpha)
            if(stem): plt.stem(y, A_plt, linefmt=stemline, markerfmt=stemmarker)
            plt.axis('tight'); plt.grid(True)
            plt.xlabel(r'$y$', fontsize=16); plt.ylabel(r'$\left|A_{n}\right|$', fontsize=16)
            if(show): plt.title(r'$\mathrm{Array}\ \mathrm{Excitation}$', fontsize=18); plt.show()
        elif (N == 1):  # i.e, linear array is along the x-direction
            plt.plot(x, A_plt, color=color, linewidth=linewidth,
                         linestyle=linestyle, alpha=alpha)
            if(stem): plt.stem(x, A_plt, linefmt=stemline, markerfmt=stemmarker)
            plt.axis('tight'); plt.grid(True)
            plt.xlabel(r'$x$', fontsize=16); plt.ylabel(r'$\left|A_{m}\right|$', fontsize=16)
            if(show): plt.title(r'$\mathrm{Array}\ \mathrm{Excitation}$', fontsize=18); plt.show()
        else:
            if (mayavi_app):  # this option opens the 3D plot in MayaVi Application
                mlab.options.backend = 'envisage'
            mlab.figure(fgcolor=fgcolor, bgcolor=bgcolor)
            s1 = mlab.quiver3d(x, y, z, z, z, A_plt)  # stem3D representation
            ranges1 = [x.min(), x.max(), y.min(), y.max(), A_plt.min(), A_plt.max()]
            mlab.axes(xlabel="x", ylabel="y", zlabel="Amn", ranges=ranges1, nb_labels=3)
            mlab.colorbar(orientation="vertical", nb_labels=5)
            s1.scene.isometric_view()
            if(show): mlab.show()
    return array_ip

def ATE(array_ip):
    r"""
    A simple function to evaluate the array taper efficiency (ATE). For the 
    definition of ATE, see [Bhattacharyya]_.

    :param array_ip: array excitation data in 'Arraytool' input format (see :func:`ip_format`)
    
    :rtype:          ATE, a Numpy float    
    """
    A = abs(array_ip[:, 3])
    A2 = A * A
    ATE = (A.sum()) ** 2 / (len(A) * A2.sum())
    return ATE

def K_norm(array_ip):
    r"""
    Function to get the normalized bore-sight slope for difference patterns. For
    the definition of K_norm, see [zinka2]_.
    
    :param array_ip: array excitation data in 'Arraytool' input format (see :func:`ip_format`)
    
    :rtype:          K_norm, a Numpy float              
    """
    X = array_ip[:, 0]
    A = array_ip[:, 3]
    AX = A * X
    A1 = abs(A)
    A2 = A1 * A1    
    K_norm = abs(AX.sum()) / np.sqrt(A2.sum())
    return K_norm

def AF_zeros(a, M, R, dist_type, nbar=False, alpha=0):
    r"""
    This function gives array-factor zeros corresponding to different
    types of array distributions. Unless you know what you are doing exactly,
    do not use this function directly. Instead, user can use the function :func:`dist`.
    
    :param a:        separation between the elements along the x-axis in wavelengths
    :param M:        number of elements along the x-axis
    :param R:        side-lobe ratio in linear scale
    :param dist_type:type of the distribution, e.g., 'Dolph-Chebyshev', 'Riblet', etc.
    :param nbar:     transition index for dilation
    :param alpha:    Taylor's asymptotic tapering parameter
    
    :rtype:          U0, a Numpy array of size (*,1)
    """
    k = 2 * np.pi  # (angular) wave-number, which is 2*pi when lambda = 1
    m = np.ceil((M - 2) / 2)
    n = np.arange(1, 1 + m, 1)  # number of zeros for symmetric array-factors
    na = np.arange(1, M, 1)  # number of zeros for 'asymmetric' array-factors
    
    if(dist_type == "Dolph-Chebyshev"):  # Dolph zeros
        c = np.cosh(np.arccosh(R) / (M - 1))
        print c
        U0 = (2 / (a * k)) * np.arccos((np.cos(np.pi * (2 * n - 1) / (2 * M - 2))) / c)        
    elif(dist_type == "Riblet"):  # Riblet zeros
        c1 = np.cosh(np.arccosh(R) / m)
        c = np.sqrt((1 + c1) / (2 + (c1 - 1) * np.cos(k * a / 2) ** 2))
        alph = c * np.cos(k * a / 2)
        xi = (1 / c) * np.sqrt(((1 + alph ** 2) / 2) + ((1 - alph ** 2) / 2) * 
                               np.cos(((2 * n - 1) * np.pi) / (2 * m)))
        U0 = (2 / (a * k)) * np.arccos(xi)
    elif(dist_type == "Duhamel-b"):  # Duhamel bi-directional end-fire array zeros
        if(a < 0.5): c = np.cosh(np.arccosh(R) / (M - 1)) / np.sin((k * a) / 2)
        else: c = np.cosh(np.arccosh(R) / (M - 1))
        U0 = (2 / (a * k)) * np.arcsin((np.cos(np.pi * (2 * n - 1) / (2 * M - 2))) / c)
    elif(dist_type == "Duhamel-u"):  # Duhamel uni-directional end-fire array zeros
        Lamb = np.cosh(np.arccosh(R) / (M - 1))
        xi = (2 / a) * (0.5 * np.pi - (np.arctan(np.tan(k * a / 2) * ((Lamb + 1) / (Lamb - 1)))))
        c = 1 / (np.sin((xi - k) * a / 2))
        U0 = -(xi / k) + (2 / (a * k)) * np.arcsin((
                            np.cos(np.pi * (2 * na - 1) / (2 * M - 2))) / c)
    elif(dist_type == "McNamara-s"):  # McNamara-Zolotarev sum-pattern zeros
        U0 = "Yet to be done"
    elif(dist_type == "McNamara-d"):  # McNamara-Zolotarev difference-pattern zeros
        if(a < 0.5): c = 1 / np.sin((k * a) / 2)
        else: c = 1
        m1 = Zol.z_m_frm_R(M - 1, R)
        xn = Zol.z_Zolotarev_poly(N=M - 1, m=m1)[1][m + 1:]
        U0 = (2 / (a * k)) * np.arcsin(xn / c)
    if(nbar):  # Taylor's Dilation procedure
        if((dist_type == "Dolph-Chebyshev") or (dist_type == "Riblet") or (dist_type == "McNamara-s")):
            n_gen = np.arange(nbar, 1 + m, 1)  # indices of the generic zeros
            U0_gen = (n_gen + alpha / 2) * (1 / (M * a))  # generic sum zeros
        elif(dist_type == "McNamara-d"):
            n_gen = np.arange(nbar, 1 + m, 1)  # indices of the generic zeros
            U0_gen = (n_gen + (alpha + 1) / 2) * (1 / (M * a))  # generic difference zeros
        sigma = U0_gen[0] / U0[nbar - 1]  # Dilation factor
        U0 = np.hstack((sigma * U0[0:nbar - 1], U0_gen))  # Dilated zeros
    U0 = np.reshape(U0, (len(U0), -1))
    return U0

def A_frm_zeros(U0, a, M, symmetry=False):
    r"""
    This function gives array excitation coefficients corresponding to the
    given array factor zeros. Unless you know what you are doing exactly, do not
    use this function directly. Instead, user can use the function :func:`dist`.

    :param U0:       arrayfactor zeros... in the format of 'AF_zeros' output form
    :param a:        separation between elements along the x-axis in wavelengths
    :param M:        number of elements along the x-axis
    :param symmetry: symmetry information to simplify numerical process... even/odd/False
    
    :rtype:          A_tot, a Numpy array of size (1,M)
    """
    k = 2 * np.pi  # (angular) wave-number, which is 2*pi when lambda = 1
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

def dist(a, M, R_x, dist_type_x, b=None, N=None, R_y=None, dist_type_y=False, mbar=False,
         nbar=False, alpha_x=0, alpha_y=0):
    r"""
    This function gives array excitation coefficients corresponding to various
    array distribution types such as Dolph-Chebyshev, McNamara-Zolotarev-sum,
    McNamara-Zolotarev-diff-f, McNamara-Zolotarev-diff-s, Taylor, Bayliss,
    Pritchard-Chebyshev-be, Pritchard-Chebyshev-ue, etc.

    :param a:           separation between the elements along the x-axis in wavelengths
    :param M:           number of elements along the x-axis
    :param R_x:         side-lobe ratio in linear scale
    :param dist_type_x: type of the distribution, e.g., 'Dolph' for Dolph-Chebyshev
    :param mbar:        transition index for dilation
    :param alpha_x:     Taylor's asymptotic tapering parameter
    
    All other parameters are similar to the above ones ... except that they
    correspond to the y-axis "principle plane" distribution.
    
    :rtype:             A_tot, a Numpy array of size (N,M)
    """
    # modify the symmetry thing in MZ-s, be, ue patterns ...
    if(dist_type_x == "Dolph-Chebyshev"):
        U0 = AF_zeros(a, M, R_x, dist_type="Dolph-Chebyshev")
        Ax = A_frm_zeros(U0, a, M, symmetry="even").T  # Done
    elif(dist_type_x == "McNamara-Zolotarev-sum"):
        U0 = AF_zeros(a, M, R_x, dist_type="Riblet")
        Ax = A_frm_zeros(U0, a, M, symmetry="even").T  # Done
    elif(dist_type_x == "McNamara-Zolotarev-diff-f"):
        U0 = AF_zeros(a, M, R_x, dist_type="McNamara-d")
        Ax = A_frm_zeros(U0, a, M, symmetry="odd").T  # Done
    elif(dist_type_x == "McNamara-Zolotarev-diff-s"):
        U0 = AF_zeros(a, M, R_x, dist_type="McNamara-d")
        Ax = A_frm_zeros(U0, a, M, symmetry="odd").T  # To be Modified later
    elif(dist_type_x == "Taylor"):
        U0 = AF_zeros(a, M, R_x, dist_type="Dolph-Chebyshev", nbar=mbar, alpha=alpha_x)
        Ax = A_frm_zeros(U0, a, M, symmetry="even").T  # Done
    elif(dist_type_x == "Bayliss"):
        U0 = AF_zeros(a, M, R_x, dist_type="McNamara-d", nbar=mbar, alpha=alpha_x)
        Ax = A_frm_zeros(U0, a, M, symmetry="odd").T  # Done
    elif(dist_type_x == "Pritchard-Chebyshev-be"):
        U0 = AF_zeros(a, M, R_x, dist_type="Duhamel-b")
        Ax = A_frm_zeros(U0, a, M, symmetry=False).T  # Done
    elif(dist_type_x == "Pritchard-Chebyshev-ue"):
        U0 = AF_zeros(a, M, R_x, dist_type="Duhamel-u")
        Ax = A_frm_zeros(U0, a, M, symmetry=False).T  # Done
        print Ax
        
    if(dist_type_y):
        # modify the symmetry thing in MZ-s, be, ue patterns ...
        if(dist_type_y == "Dolph-Chebyshev"):
            V0 = AF_zeros(b, N, R_y, dist_type="Dolph-Chebyshev")
            Ay = A_frm_zeros(V0, b, N, symmetry="even")  # Done
        elif(dist_type_y == "McNamara-Zolotarev-sum"):
            V0 = AF_zeros(b, N, R_y, dist_type="Riblet")
            Ay = A_frm_zeros(V0, b, N, symmetry="even")  # Done
        elif(dist_type_y == "McNamara-Zolotarev-diff-f"):
            V0 = AF_zeros(b, N, R_y, dist_type="McNamara-d")
            Ay = A_frm_zeros(V0, b, N, symmetry="odd")  # Done
        elif(dist_type_y == "McNamara-Zolotarev-diff-s"):
            V0 = AF_zeros(b, N, R_y, dist_type="McNamara-d")
            Ay = A_frm_zeros(V0, b, N, symmetry="odd")  # To be Modified later
        elif(dist_type_y == "Taylor"):
            V0 = AF_zeros(b, N, R_y, dist_type="Dolph-Chebyshev", nbar=nbar, alpha=alpha_y)
            Ay = A_frm_zeros(V0, b, N, symmetry="even")  # Done
        elif(dist_type_y == "Bayliss"):
            V0 = AF_zeros(b, N, R_y, dist_type="McNamara-d", nbar=nbar, alpha=alpha_y)
            Ay = A_frm_zeros(V0, b, N, symmetry="odd")  # Done
        elif(dist_type_y == "Pritchard-Chebyshev-be"):
            V0 = AF_zeros(b, N, R_y, dist_type="Duhamel-b")
            Ay = A_frm_zeros(V0, b, N, symmetry=False)  # Done
        elif(dist_type_y == "Pritchard-Chebyshev-ue"):
            V0 = AF_zeros(b, N, R_y, dist_type="Duhamel-u")
            Ay = A_frm_zeros(V0, b, N, symmetry=False)  # Done
        Ax = np.tile(Ax, (N, 1))
        Ay = np.tile(Ay, (1, M))
        Ax = Ax * Ay
        
    A_tot = Ax
    return A_tot

def cutoff(F, dB_limit=-40):
    r"""
    When AF/GF/NF is 0, their dB value is '-infinity'. So, this function
    will be used to cut-off all the value below some 'dB_limit'.
    
    :param F:            F is a Numpy array (in our case, it is usually AF/GF/NF)
    :param dB_limit:     cut-off level in dB, default value is -40
    
    :rtype:              a Numpy array, same size as the input parameter F
    """
    msk1 = F < dB_limit
    fill = msk1 * dB_limit
    msk2 = F >= dB_limit
    F = F * (msk2) + fill
    return F

def pattern_u(array_ip, u_scan=0, u_min=-1, u_max=1, u_num=200, scale="dB",
              dB_limit=-40, factor="GF", plot_type="rect", lattice=False,
              color='b', linewidth=1, linestyle='-', alpha=1, show=True):
    r"""
    Function to evaluate 2D AF/GF/NF/NF0 of a linear array in u-domain (so, array
    should be along the x-axis). By default, this function calculates the gain-factor(GF).
    For the definitions of GF, NF, AF, go through the Arraytool documentation.
   
    :param array_ip:      array excitation data in 'Arraytool' input format (see :func:`ip_format`)
    :param u_scan:        beam scan position
    :param u_min, u_max:  limits of u-domain
    :param u_num:         number of points between 'u_min' and 'u_max' including 
                          boundaries
    :param scale:         specifies the scale choice ... dB/linear
    :param dB_limit:      cutoff limit (see :func:`cutoff`)
    :param factor:        type of pattern you need ... AF/GF/NF/NF0
    :param plot_type:     can be rect/polar ... if False, nothing happens
    :param lattice:       If True, highlights visible-space and lattice period in 
                          "rect" plot mode
    
    Lattice period is meaningful only if the array is "uniformly spaced". 
    
    All other parameters are nothing but 'Matplotlib' parameters. These should 
    be familiar to 'Matlab' or 'Matplotlib' users.
    
    :rtype:               [u,F] ... i.e., both u and F are column matrices ...
                          F can be AF/GF/NF/NF0 ... and can be either in linear
                          or in dB scale ...
    """
    x = array_ip[:, 0]
    y = array_ip[:, 1]
    z = array_ip[:, 2]
    A = array_ip[:, 3]  # un-packing "array_ip" finished
    k = 2 * np.pi  # (angular) wave-number, which is 2*pi when lambda = 1

    # Making sure all elements in y and z columns of the "array_ip" are zeros
    z_flag = True # z_flag, i.e., zinka's flag ... :)
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

        # Evaluating array-factor of the linear array-2
        AF = np.dot(np.exp(1j * k * U * X), A)

        # Evaluation of F = (AF/GF/NF) => depending upon the user's choice
        if(factor == "AF"):
            F = AF; n1 = ""; ff = r"$\mathrm{Array-Factor}\ $"; f1 = r"$AF\ $"
        elif(factor == "GF"):
            P_inc = ((abs(A)) ** 2).sum()
            GF = AF / np.sqrt(P_inc)  # Converting the AF to GF ... AKB's definition
            F = GF; n1 = ""; ff = r"$\mathrm{Gain-Factor}\ $"; f1 = r"$GF\ $"
        elif(factor == "NF"):
            norm_fact = (abs(A)).sum()
            F = AF / norm_fact
            n1 = r"$\mathrm{Normalized}\ $"; ff = r"$\mathrm{Factor}\ $"; f1 = r"$NF\ $"
        elif(factor == "NF0"):
            norm_fact = AF.max()            
            F = AF / norm_fact
            n1 = r"$\mathrm{Normalized}\ $"; ff = r"$\mathrm{Factor}\ $"; f1 = r"$NF0\ $"            

        # converting 'F' from linear to dB scale, if needed
        if(scale == "linear"):
            F_plt = abs(F)
#            F_plt = F # temporarily ... delete this later ... can be problematic in complex array factor cases
            ss = r"$\mathrm{in}\ \mathrm{linear}\ \mathrm{scale}$"
        elif(scale == "dB"):
            F = 20 * np.log10(abs(F))
            # cutoff the "F" below some limit ... just for the plotting purpose
            F_plt = cutoff(F, dB_limit)
            ss = r"$\mathrm{in}\ \mathrm{dB}\ \mathrm{scale}$"

        # plotting the factor 'F_plt'
        if(plot_type):
            if(plot_type == "rect"):  # rectangular plot
                plt.plot(u, F_plt, color=color, linewidth=linewidth,
                         linestyle=linestyle, alpha=alpha)
                if(lattice):  # highlighting the visible-space and unit-lattice
                    plt.axvspan(-1, +1, facecolor='y', alpha=0.2)
                    lim = -np.pi / ((x[2] - x[1]) * k)
                    plt.axvspan(-lim, +lim, facecolor='b', alpha=0.2)
                plt.axis('tight'); plt.grid(True)
                plt.xlabel(r'$u,\ \mathrm{where}\ u=\sin \theta\ \mathrm{in}\ \mathrm{the}\ \mathrm{visible-space}$', fontsize=16)
                plt.ylabel(f1 + r'$(u)$', fontsize=16)
            if(plot_type == "polar"):  # polar plot
                tht = np.arcsin(u)
                
                # the below code is to make sure that tht=0 starts in the 'North'
                testFig = plt.figure('Plot')
                rect = [0.1, 0.1, 0.8, 0.8]
                testAx = testFig.add_axes(rect, projection='polar')
                testAx.set_theta_zero_location('N')
                testAx.set_theta_direction(-1)
                
                if(scale == "linear"):
                    plt.polar(tht, F_plt, color=color, linewidth=linewidth,
                              linestyle=linestyle, alpha=alpha)
                    plt.polar(np.pi - tht, F_plt, color=color, linewidth=linewidth,
                              linestyle=linestyle, alpha=alpha)
                if(scale == "dB"):
                    plt.polar(tht, F_plt - dB_limit, color=color,
                              linewidth=linewidth, linestyle=linestyle, alpha=alpha)
                    plt.polar(np.pi - tht, F_plt - dB_limit, color=color,
                              linewidth=linewidth, linestyle=linestyle, alpha=alpha)
            plt.title(n1 + ff + ss + '$\ (\phi=0^\circ)$', fontsize=18)          
            
            if(show): plt.show()
    return u, F

def pattern_uv(array_ip, u_scan=0, v_scan=0, u_min=-1, u_max=1, u_num=50,
               v_min=-1, v_max=1, v_num=50, scale="dB", dB_limit=-40, 
               factor="GF", plot_type="rect", fgcolor=(1, 1, 1),
               bgcolor=(0.5, 0.5, 0.5), mayavi_app=False):
    r"""
    Function to evaluate 3D AF/GF/NF/NF0 of a planar array in uv-domain.
    By default, this function calculates the gain-factor (GF).
   
    :param array_ip:       array excitation data in 'Arraytool' input format (see :func:`ip_format`)
    :param u_scan, v_scan: beam scan position in uv-domain
    :param u_min, etc:     limits of uv-domain
    :param u_num, v_num:   number of points between 'u_min' and 'u_max' including boundaries
    :param scale:          specifies the scale choice ... dB/linear
    :param dB_limit:       cutoff limit (see :func:`cutoff`)
    :param factor:         type of pattern you need ... AF/GF/NF/NF0
    :param plot_type:      can be rect/contour ... if False, nothing happens
    :param mayavi_app:     if True, the 3D plot will be opened in the MayaVi application
    
    :rtype:                [u,v,F] ... F can be AF/GF/NF/NF0 ... 
                           and can be either in linear or in dB scale ...
    """
    x = array_ip[:, 0]
    y = array_ip[:, 1]
    z = array_ip[:, 2]
    A = array_ip[:, 3]  # un-packing "array_ip" finished
    k = 2 * np.pi  # (angular) wave-number, which is 2*pi when lambda = 1
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
            GF = AF / np.sqrt(P_inc)  # Converting the AF to GF
            F = GF; n1 = ""; ff = "Gain-Factor "; f1 = "GF "
        elif(factor == "NF"):
            norm_fact = (abs(A)).sum()
            F = AF / norm_fact
            n1 = "Normalized "; ff = "Factor "; f1 = "NF "
        elif(factor == "NF0"):
            norm_fact = AF.max() 
            F = AF / norm_fact
            n1 = "Normalized "; ff = "Factor "; f1 = "NF0 "            

        # converting 'F' from linear to dB scale, if needed
        if(scale == "linear"):
            F_plt = abs(F)
#            F_plt = F # temporarily ... delete this later ... can be problematic in complex array factor cases            
            ss = "in linear scale"
        elif(scale == "dB"):
            F = 20 * np.log10(abs(F))
            # cutoff the "F" below some limit ... just for the plotting purpose
            F_plt = cutoff(F, dB_limit)
            ss = "in dB scale"

        # plotting F_plt
        if(plot_type):
            if(plot_type == "rect"):  # rectangular plot
                if (mayavi_app):  # opens the 3D plot in MayaVi Application
                    mlab.options.backend = 'envisage'                
                mlab.figure(fgcolor=fgcolor, bgcolor=bgcolor)
                plt3d = mlab.surf(u, v, F_plt, warp_scale='auto')
                ranges1 = [u_min, u_max, v_min, v_max, F_plt.min(), F_plt.max()]
                mlab.axes(xlabel='u', ylabel='v', zlabel=f1,
                          ranges=ranges1, nb_labels=5)
                mlab.title(n1 + ff + ss, size=0.35)
                mlab.colorbar(orientation="vertical", nb_labels=5)
                plt3d.scene.isometric_view()
                mlab.show()                
            if(plot_type == "contour"):  # contour plot
                plt.contourf(u, v, F_plt)
                vs = plt.Circle((0, 0), radius=1, edgecolor='w', fill=False)
                ax = plt.gca(); ax.add_patch(vs)                
                plt.axis('image'); plt.grid(True)
                plt.xlabel(r'$u,\ \mathrm{where}\ u=\sin \theta \cos \phi\ \mathrm{in}\ \mathrm{the}\ \mathrm{visible-space}$', fontsize=16)
                plt.ylabel(r'$v,\ \mathrm{where}\ v=\sin \theta \sin \phi\ \mathrm{in}\ \mathrm{the}\ \mathrm{visible-space}$', fontsize=16)
                plt.colorbar(format='$%.2f$')
                plt.show()                
    return u, v, F

def pattern_tp(array_ip, tht_scan=0, phi_scan=0, tht_min=0, tht_max=np.pi, tht_num=50,
               phi_min=0, phi_max=2 * np.pi, phi_num=50, scale="dB",
               dB_limit=-40, factor="GF", plot_type="polar", fgcolor=(1, 1, 1),
               bgcolor=(0.5, 0.5, 0.5), mayavi_app=False):
    r"""
    Function to evaluate 3D AF/GF/NF/NF0 of an arbitrary 3D array in (tht, phi)-domain.
    By default, this function calculates the gain-factor (GF).
   
    :param array_ip:       array excitation data in 'Arraytool' input format (see :func:`ip_format`)
    :param tht_scan, etc:  beam scan position in (tht, phi)-domain
    :param tht_min, etc:   limits of (tht, phi)-domain
    :param tht_num, etc:   number of points between 'tht_min' and 'tht_max' including
                           the boundaries
    :param scale:          specifies the scale choice ... dB/linear
    :param dB_limit:       cutoff limit (see :func:`cutoff`)
    :param factor:         type of pattern you need ... AF/GF/NF/NF0
    :param plot_type:      can be rect/polar/contour ... if False, nothing happens
    :param mayavi_app:     if True, the 3D plot will be opened in the MayaVi application
    
    :rtype:                [tht,phi,F] ... F can be AF/GF/NF/NF0 ... 
                           and can be either in linear or in dB scale ...
    """
    x = array_ip[:, 0]
    y = array_ip[:, 1]
    z = array_ip[:, 2]
    A = array_ip[:, 3]  # un-packing "array_ip" finished
    k = 2 * np.pi  # (angular) wave-number, which is 2*pi when lambda = 1
    tht_numj = complex(0, tht_num)
    phi_numj = complex(0, phi_num)

    [tht, phi] = np.mgrid[tht_min:tht_max:tht_numj, phi_min:phi_max:phi_numj]
    u = np.sin(tht) * np.cos(phi); v = np.sin(tht) * np.sin(phi); w = np.cos(tht)
    u1 = np.reshape(u, (u.size, -1))
    v1 = np.reshape(v, (v.size, -1))
    w1 = np.reshape(w, (w.size, -1))
    u_scan = np.sin(tht_scan) * np.cos(phi_scan)
    v_scan = np.sin(tht_scan) * np.sin(phi_scan)
    w_scan = np.cos(tht_scan)
#    w_scan = 0 #delete this later
    
    A = np.reshape(A, (len(A), -1))
    U = np.tile(u1 - u_scan, len(x))
    V = np.tile(v1 - v_scan, len(x))
    W = np.tile(w1 - w_scan, len(x))
    X = np.tile(x, (u.size, 1))
    Y = np.tile(y, (u.size, 1))
    Z = np.tile(z, (u.size, 1))

    # Evaluating array-factor of the planar array
    AF1 = np.dot(np.exp(1j * k * (U * X + V * Y + W * Z)), A)
    AF = np.reshape(AF1, u.shape)

    # Evaluation of F = (AF/GF/NF) => depending upon the user's choice
    if(factor == "AF"):
        F = AF; n1 = ""; ff = "Array-Factor "; f1 = "AF "
    elif(factor == "GF"):
        P_inc = ((abs(A)) ** 2).sum()
        GF = AF / np.sqrt(P_inc)  # Converting the AF to GF
        F = GF; n1 = ""; ff = "Gain-Factor "; f1 = "GF "
    elif(factor == "NF"):
        norm_fact = (abs(A)).sum()
        F = AF / norm_fact
        n1 = "Normalized "; ff = "Factor "; f1 = "NF "
    elif(factor == "NF0"):
        norm_fact = AF.max() 
        F = AF / norm_fact
        n1 = "Normalized "; ff = "Factor "; f1 = "NF0 "         

    # converting 'F' from linear to dB scale, if needed
    if(scale == "linear"):
        F_plt = abs(F)
#         F_plt = F # temporarily ... delete this later ... can be problematic in complex array factor cases                    
        ss = "in linear scale"
    elif(scale == "dB"):
        F = 20 * np.log10(abs(F))
        # cutoff the "F" below some limit ... just for the plotting purpose
        F_plt = cutoff(F, dB_limit)
        ss = "in dB scale"

    # plotting the factor (AF/GF/NF)
    if(plot_type):
        if (mayavi_app):  # opens the 3D plot in MayaVi Application
            mlab.options.backend = 'envisage'
        if(plot_type == "rect"):  # rectangular plot
            mlab.figure(fgcolor=fgcolor, bgcolor=bgcolor)
            plt3d = mlab.surf(tht, phi, F_plt, warp_scale='auto')
            ranges1 = [tht.min(), tht.max(), phi.min(), phi.max(), F_plt.min(), F_plt.max()]
            mlab.axes(xlabel='Tht', ylabel='Phi', zlabel=f1,
                      ranges=ranges1, nb_labels=5)
            mlab.title(n1 + ff + ss, size=0.35)
            mlab.colorbar(orientation="vertical", nb_labels=5)
            plt3d.scene.isometric_view()
            mlab.show()
        if(plot_type == "polar"):  # rectangular plot
            if(scale == "dB"):
                F_plt = F_plt - dB_limit
            F_plt_x = F_plt * u; F_plt_y = F_plt * v; F_plt_z = F_plt * w            
            ranges1 = [F_plt_x.min(), F_plt_x.max(), F_plt_y.min(), F_plt_y.max(), F_plt_z.min(), F_plt_z.max()]
            mlab.figure(fgcolor=fgcolor, bgcolor=bgcolor)
            plt3d = mlab.mesh(F_plt_x, F_plt_y, F_plt_z, scalars=F_plt, extent=ranges1)
            mlab.axes(xlabel='x', ylabel='y', zlabel='z',
                      ranges=ranges1, nb_labels=5)
            mlab.title(n1 + ff + ss, size=0.35)
            mlab.colorbar(orientation="vertical", nb_labels=5)
            plt3d.scene.isometric_view()
            mlab.show()        
        if(plot_type == "contour"):  # contour plot
            plt.contourf(tht, phi, F_plt)
            plt.axis('tight'); plt.grid(True)
            plt.xlabel(r'$\theta$', fontsize=16)
            plt.ylabel(r'$\phi$', fontsize=16)
            plt.colorbar(format='$%.2f$')
            plt.show()                
    return tht, phi, F

def pattern_t(array_ip, tht_scan=0, phi_scan=0, tht=0.5 * np.pi,
               phi_min=0, phi_max=2 * np.pi, phi_num=200, scale="dB",
               dB_limit=-40, factor="GF", plot_type="rect", color='b',
               linewidth=1, linestyle='-', alpha=1, show=True):
    r"""
    Function to evaluate 2D AF/GF/NF/NF0 of an arbitrary 3D array (theta = constant).
    By default, this function calculates the gain-factor (GF).
   
    :param array_ip:       array excitation data in 'Arraytool' input format (see :func:`ip_format`)
    :param tht_scan, etc:  beam scan position in (tht, phi)-domain
    :param tht:            'theta = constant' cut specification
    :param phi_min, etc:   limits of phi-domain
    :param phi_num:        number of points between 'phi_min' and 'phi_max' 
                           including the boundaries
    :param scale:          specifies the scale choice ... dB/linear
    :param dB_limit:       cutoff limit (see :func:`cutoff`)
    :param factor:         type of pattern you need ... AF/GF/NF/NF0
    :param plot_type:      can be rect/polar ... if False, nothing happens
    
    All other parameters are nothing but 'Matplotlib' parameters. These should 
    be familiar to 'Matlab' or 'Matplotlib' users.    
    
    :rtype:                [phi,F] ... F can be AF/GF/NF/NF0 ... 
                           and can be either in linear or in dB scale ...
    """    
    tht, phi, F = pattern_tp(array_ip, tht_scan=tht_scan, phi_scan=phi_scan, tht_min=tht, tht_max=tht, tht_num=1,
               phi_min=phi_min, phi_max=phi_max, phi_num=phi_num, scale=scale,
               dB_limit=dB_limit, factor=factor, plot_type=False)
    
    # converting 'F' to absolute value or cutting off the dB values
    if(scale == "linear"):
        F_plt = abs(F) # pattern_tp fives linear values ... some of them can be -ve
    elif(scale == "dB"):
        F_plt = cutoff(F, dB_limit) # pattern_tp gives F in dB without any cutoff
    
    # For title
    if(factor == "AF"):
        ff = r"$\mathrm{Array-Factor}\ $"
    elif(factor == "GF"):
        ff = r"$\mathrm{Gain-Factor}\ $"
    elif(factor == "NF"):
        ff = r"$\mathrm{Normalized-Factor}\ NF\  $"
    elif(factor == "NF0"):
        ff = r"$\mathrm{Normalized-Factor}\ (NF0)\ $"        
    if(scale == "linear"):
        ss = r"$\mathrm{in}\ \mathrm{linear}\ \mathrm{scale}$"
    elif(scale == "dB"):
        ss = r"$\mathrm{in}\ \mathrm{dB}\ \mathrm{scale}$"    
    
    # plotting the factor 'F_plt'
    if(plot_type):
        if(plot_type == "rect"):  # rectangular plot
            plt.plot(phi[0], F_plt[0], color=color, linewidth=linewidth,
                     linestyle=linestyle, alpha=alpha)
            plt.axis('tight'); plt.grid(True)
            plt.xlabel(r'$\phi$', fontsize=16)
            plt.ylabel(r'$AF$', fontsize=16)
        if(plot_type == "polar"):  # polar plot
            if(scale == "linear"):
                plt.polar(phi[0], F_plt[0], color=color, linewidth=linewidth,
                          linestyle=linestyle, alpha=alpha)
            if(scale == "dB"):
                plt.polar(phi[0], F_plt[0] - dB_limit, color=color,
                          linewidth=linewidth, linestyle=linestyle, alpha=alpha)
        plt.title(ff + ss, fontsize=18)
        if(show): plt.show()
    return phi[0], F[0]

def pattern_p(array_ip, tht_scan=0, phi_scan=0, phi=0 * np.pi,
               tht_min=0, tht_max=1 * np.pi, tht_num=200, scale="dB",
               dB_limit=-40, factor="GF", plot_type="rect", color='b',
               linewidth=1, linestyle='-', alpha=1, show=True):
    r"""
    Function to evaluate 2D AF/GF/NF of an arbitrary 3D array (phi = constant).
    By default, this function calculates the gain-factor (GF).
   
    :param array_ip:       array excitation data in 'Arraytool' input format (see :func:`ip_format`)
    :param tht_scan, etc:  beam scan position in (tht, phi)-domain
    :param phi:            'phi = constant' cut specification
    :param tht_min, etc:   limits of tht-domain
    :param tht_num:        number of points between 'tht_min' and 'tht_max' 
                           including the boundaries
    :param scale:          specifies the scale choice ... dB/linear
    :param dB_limit:       cutoff limit (see :func:`cutoff`)
    :param factor:         type of pattern you need ... AF/NF/GF
    :param plot_type:      can be rect/polar ... if False, nothing happens
    
    All other parameters are nothing but 'Matplotlib' parameters. These should 
    be familiar to 'Matlab' or 'Matplotlib' users.    
    
    :rtype:                [tht,F] ... F can be AF/GF/NF/NF0 ... 
                           and can be either in linear or in dB scale ...
    """
    tht, phi, F = pattern_tp(array_ip, tht_scan=tht_scan, phi_scan=phi_scan,
                             tht_min=tht_min, tht_max=tht_max, tht_num=tht_num,
                             phi_min=phi, phi_max=phi, phi_num=1, scale=scale,
                             dB_limit=dB_limit, factor=factor, plot_type=False)

    # converting 'F' to absolute value or cutting off the dB values
    if(scale == "linear"):
        F_plt = abs(F) # pattern_tp fives linear values ... some of them can be -ve
    elif(scale == "dB"):
        F_plt = cutoff(F, dB_limit) # pattern_tp gives F in dB without any cutoff    
    
    # For title
    if(factor == "AF"):
        ff = r"$\mathrm{Array-Factor}\ $"
    elif(factor == "GF"):
        ff = r"$\mathrm{Gain-Factor}\ $"
    elif(factor == "NF"):
        ff = r"$\mathrm{Normalized-Factor}\ NF\ $"
    elif(factor == "NF0"):
        ff = r"$\mathrm{Normalized-Factor}\ (NF0)\ $"
    if(scale == "linear"):
        ss = r"$\mathrm{in}\ \mathrm{linear}\ \mathrm{scale}$"
    elif(scale == "dB"):
        ss = r"$\mathrm{in}\ \mathrm{dB}\ \mathrm{scale}$"    
    
    # plotting the factor 'F_plt'
    if(plot_type):
        if(plot_type == "rect"):  # rectangular plot
            plt.plot(tht[:, 0], F[:, 0], color=color, linewidth=linewidth,
                     linestyle=linestyle, alpha=alpha)
            plt.axis('tight'); plt.grid(True)
            plt.xlabel(r'$\phi$', fontsize=16)
            plt.ylabel(r'$AF$', fontsize=16)
        if(plot_type == "polar"):  # polar plot
            
            # the below code is to make sure that tht=0 starts in the 'North'
            testFig = plt.figure('Plot')
            rect = [0.1, 0.1, 0.8, 0.8]
            testAx = testFig.add_axes(rect, projection='polar')
            testAx.set_theta_zero_location('N')
            testAx.set_theta_direction(-1)            
            
            if(scale == "linear"):
                plt.polar(tht[:, 0], F[:, 0], color=color, linewidth=linewidth,
                          linestyle=linestyle, alpha=alpha)
            if(scale == "dB"):
                plt.polar(tht[:, 0], F[:, 0] - dB_limit, color=color,
                          linewidth=linewidth, linestyle=linestyle, alpha=alpha)
        plt.title(ff + ss, fontsize=18)
        if(show): plt.show()
    return tht[:, 0], F[:, 0]

if __name__ == '__main__':

    # frequency and array-arrangement (actual values)
    freq = 10e9  # frequency of operation in Hzs
    wav_len = 3e8 / freq  # wavelength in meters
    M = 10  # no. of elements along the x-axis
    N = 1  # no. of elements along the y-axis
    a1 = 15e-3  # separation between elements along the x-axis in meters
    b1 = 17e-3  # separation between elements along the y-axis in meters
    gamma = np.pi / 2  # lattice angle in radians

    # normalized values
    a = a1 / wav_len  # 'a1' in-terms of lambda (wavelength)
    b = b1 / wav_len  # 'b1' in-terms of lambda (wavelength)

    SLR = 25  # side-lobe ratio in dB
    R = 10 ** (SLR / 20)  # converting SLR from dB scale to linear scale

#     A = np.array([[1,2,3,4,5]]) # manually entering the coefficients
    A = np.ones((N, M))  # Uniform excitation     
#    A = np.random.rand(N, M) # Random excitation
 
#    # Using the function 'AF_zeros' to find arrayfactor zeros
#    U0 = AF_zeros(a, M, R, dist_type="Dolph", nbar=5, alpha=0)
#    print 'arrayfactor zeros:', '\n', U0
 
#    # Obtaining array excitation coefficients from the arrayfactor zeros
#    A = A_frm_zeros(U0, a, M, symmetry="even").T # finding excitation coefficients
#    print 'array coefficients:', '\n', A.T
 
    # Finding the array excitation directly from "dist" function
#     A = dist(a, M, R, dist_type_x="McNamara-Zolotarev-diff-f", mbar=False, alpha_x=0)
#     print 'array coefficients:', '\n', A.T
 
#    # Converting the 'excitation & position' info into 'Arraytool' input format
    array_ip = ip_format(a, b, A, gamma, plot=False, stem=True, mayavi_app=False)
    
#    ATE = ATE(array_ip) # array taper efficiency;
#    print type(ATE)
 
# #     Calling the 'pattern_u' function to evaluate and plot 2D AF/GF/NF/NF0
#     u,F = pattern_u(array_ip, u_scan=0, u_min=-1, u_max=1, u_num=700, scale="dB",
#               dB_limit=-40, factor="NF0", plot_type="polar", lattice=True)
 
#    # Calling the 'pattern_uv' function to evaluate and plot 3D AF/GF/NF/NF0
#     pattern_uv(array_ip, u_scan=0, v_scan=0, u_min= -3, u_max=3, u_num=100,
#                v_min= -3, v_max=3, v_num=100, scale="dB", dB_limit=-40,
#                factor="NF0", plot_type="rect", mayavi_app=False)
     
#     Calling the 'pattern_tp' function to evaluate and plot 3D AF/GF/NF/NF0  
#     pattern_tp(array_ip, tht_scan=(0)*np.pi, phi_scan=(0)*np.pi, tht_min= 0, tht_max=1*np.pi, tht_num=100,
#                phi_min= 0*np.pi, phi_max=2*np.pi, phi_num=100, scale="dB", dB_limit= -40,
#                factor="NF0", plot_type="polar", mayavi_app=False)

#     Calling the 'pattern_t' function to evaluate and plot 2D AF/GF/NF/NF0 
#     pattern_t(array_ip, tht_scan=0, phi_scan=0, tht=0.5 * np.pi,
#                    phi_min=0, phi_max=2 * np.pi, phi_num=200, scale="dB",
#                    dB_limit=-40, factor="NF0", plot_type="polar", color='b',
#                    linewidth=1, linestyle='-', alpha=1, show=True)
    
#     Calling the 'pattern_p' function to evaluate and plot 2D AF/GF/NF/NF0     
    pattern_p(array_ip, tht_scan=0, phi_scan=0, phi=0 * np.pi,
                   tht_min=0, tht_max=2 * np.pi, tht_num=200, scale="dB",
                   dB_limit=-40, factor="NF0", plot_type="polar", color='b',
                   linewidth=1, linestyle='-', alpha=1, show=True)    
 
#==============================================================================
# Programming tasks (NOTES to myself)
#==============================================================================

# odd symmetry ... optimize A_from_zeros
# also A_from_zeros and dist are not fully done ... implement a few other important patterns
# use same names in AF_zeros and dist

# use backslah instead of inverse
