#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2014 Srinivasa Rao Zinka
# License: New BSD License.

from __future__ import division
import numpy as np
from scipy import integrate, special
import matplotlib.pyplot as plt
import planar as pl
import Zolotarev as zl
from mayavi import mlab
import warnings

# adjusting "matplotlib" label fonts
from matplotlib import rc
rc('text', usetex=True)

def ip_format_c(N, radius, A="uniform", starts_at_zero=True, plot_type="2D",
                color='b', linewidth=1, linestyle='-', alpha=1, show=True,
                stem=False, stemline='g--', stemmarker='ro', fgcolor=(1, 1, 1),
                bgcolor=(0.5, 0.5, 0.5), mayavi_app=False):
    r"""
    Function to generate the 'Arraytool' input format for circular ring arrays.

    :param N:          number of elements in the uniformly spaced circular ring array
    :param radius:     radius of the circular ring in wavelengths
    :param A:          a 'column matrix' specifying the excitation values of the 
                       circular ring array; by default it will be uniform excitation
    :param plot_type:  can be '2D'/'3D' ... if False, nothing happens    
    :param stem:       if True, the array excitation is plotted as 'stem plot'
    :param mayavi_app: if True, the 3D plot will be opened in the MayaVi application
    :param starts_at_zero:   'True' if array starts at beta=0    
    
    All other parameters are nothing but the 'Matplotlib/Mayavi' parameters. 
    These should be familiar to 'Matlab' or 'Matplotlib/Mayavi' users.
    
    :rtype:            array_ip, a Numpy array of size (Number of elements(A),4)
    """
    
    # Creating Arraytool input form 'array_ip' for the circular ring array
    if (A == "uniform"):
        A = np.ones((N, 1))
    if (starts_at_zero):
        position_beta = (np.linspace(1, N, num=N) - 1) * (2 * np.pi / N)
    else:
        position_beta = (np.linspace(1, N, num=N) - 0.5) * (2 * np.pi / N)
    position_beta = np.reshape(position_beta, (N, -1))
    position_x = radius * np.cos(position_beta)
    position_y = radius * np.sin(position_beta)
    position_z = np.zeros_like(position_x)
    array_ip = np.hstack((position_x, position_y, position_z, A))
    
    # Plotting 2D/3D plots
    if (plot_type):
        # checking whether 'A' has any imaginary values
        if((abs(A.imag) > 1e-10).sum()):
            A_plt = abs(A)  # if A.imag are significant, then '|A|' will be plotted
            warnings.warn('Since, the given excitation "A" has significant imaginary parts, stem plot for abs(A) is plotted')
        else:
            A_plt = A.real  # if A.imag are negligible, then 'A'  will be plotted
            warnings.warn('Since, the given excitation "A" has very small imaginary parts, stem plot for "A.real" is plotted')
        if (plot_type == "2D"):  # plot 2D plot in Matplotlib
            plt.plot(position_beta, A_plt, color=color, linewidth=linewidth,
                         linestyle=linestyle, alpha=alpha)
            if(stem): plt.stem(position_beta, A_plt, linefmt=stemline, markerfmt=stemmarker)
            plt.axis('tight'); plt.grid(True)
            plt.xlabel(r'$y$', fontsize=16); plt.ylabel(r'$\left|A_{n}\right|$', fontsize=16)
            if(show): plt.title(r'$\mathrm{Array}\ \mathrm{Excitation}$', fontsize=18); plt.show()
        else:
            if (mayavi_app):  # this option opens the 3D plot in MayaVi Application
                mlab.options.backend = 'envisage'
            mlab.figure(fgcolor=fgcolor, bgcolor=bgcolor)
            s1 = mlab.quiver3d(position_x, position_y, position_z, position_z, position_z, A_plt)  # stem3D representation
            ranges1 = [position_x.min(), position_x.max(), position_y.min(), position_y.max(), A_plt.min(), A_plt.max()]
            mlab.axes(xlabel="x", ylabel="y", zlabel="A", ranges=ranges1, nb_labels=3)
            mlab.colorbar(orientation="vertical", nb_labels=5)
            s1.scene.isometric_view()
            if(show): mlab.show()

    return array_ip

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

def eval_Taylor(P, R, mbar, alpha_x, x):
    """My function description here.""" 
    
    if(P%2==0):
        T0 = 2*np.pi
    else:
        T0=1*np.pi
    
    A = pl.dist(1, P + 1, R, dist_type_x='Taylor', mbar=mbar, alpha_x=alpha_x)
    x, result = IFS(A, T0, m_start= -P / 2, m_stop=P / 2, x_min=x, x_max=x, x_num=1)
            
#    if(P%2==0):
#        A = pl.dist(1, P + 1, R, dist_type_x='Taylor', mbar=mbar, alpha_x=alpha_x)
#        x, result = IFS(A, T0=2*np.pi, m_start= -P / 2, m_stop=P / 2, x_min=x, x_max=x, x_num=1)
#    else:
#        A = pl.dist(1, 2*P + 1, R, dist_type_x='Taylor', mbar=mbar, alpha_x=alpha_x)
#        x, result = IFS(A, T0=1*np.pi, m_start= -P, m_stop=P, x_min=x, x_max=x, x_num=1)
    
    result = result[0,0].real
    
    return result

def eval_Bayliss(P, R, mbar, alpha_x, x):
    """My function description here.""" 
    
    if(P%2==0):
        print "Order needs to be an ODD number for null patterns"
    else:
        T0=1*np.pi
    
    A = pl.dist(1, P + 1, R, dist_type_x='Bayliss', mbar=mbar, alpha_x=alpha_x)
    x, result = IFS(A, T0, m_start= -P / 2, m_stop=P / 2, x_min=x, x_max=x, x_num=1)
    
    result = result[0,0].imag
    
    return result

def FS_Taylor(N, R, mbar, alpha_x, x_min, x_max, x_num, plot_far=False, dB_limit= -40):
    """Function to evaluate Fourier series coefficients of Chebyshev far-field
       pattern"""
       
    R = str(R)
    mbar = str(mbar)
    alpha_x = str(alpha_x)
           
    if(N % 2 == 0):
        m_start = int(-N / 2)
        m_stop = int(N / 2)
        N = str(N)
        fun_str_re = 'eval_Taylor(' + N + ',' + R + ',' + mbar + ',' + alpha_x + ',' + 'x' +')'
        print fun_str_re
        m_index, zm = FS(fun_str_re, m_start=m_start, m_stop=m_stop, err_lim=1e-5)
    else:
        m_start = -N
        m_stop = N
        N = str(N)
        fun_str_re = 'eval_Taylor(' + N + ',' + R + ',' + mbar + ',' + alpha_x + ',' + 'x' +')'
        print fun_str_re
        m_index, zm = FS(fun_str_re, m_start=m_start, m_stop=m_stop, err_lim=1e-5)
        
    if(plot_far):
        x, AF = IFS(zm, 2 * np.pi, m_start, m_stop, x_min, x_max, x_num)        
        AF = 20 * np.log10(abs(AF))
        AF = pl.cutoff(AF, dB_limit)
        plt.plot(x * (180 / np.pi), AF); plt.axis('tight'); plt.grid(True)        
        plt.title('Far-field Pattern')
        plt.xlabel(r'$\phi$')
        plt.ylabel('AF')
        plt.show()        
                     
    return m_index, zm

def FS_Bayliss(N, R, mbar, alpha_x, x_min, x_max, x_num, plot_far=False, dB_limit= -40):
    """Function to evaluate Fourier series coefficients of Chebyshev far-field
       pattern"""
       
    R = str(R)
    mbar = str(mbar)
    alpha_x = str(alpha_x)
           
    if(N % 2 == 0):
        print "Order needs to be an ODD number for null patterns"
    else:
        m_start = -N
        m_stop = N
        N = str(N)
        fun_str_re = 'eval_Bayliss(' + N + ',' + R + ',' + mbar + ',' + alpha_x + ',' + 'x' +')'
        print fun_str_re
        m_index, zm = FS(fun_str_re, m_start=m_start, m_stop=m_stop, err_lim=1e-5)
        
    if(plot_far):
        x, AF = IFS(zm, 2 * np.pi, m_start, m_stop, x_min, x_max, x_num)        
        AF = 20 * np.log10(abs(AF))
        AF = pl.cutoff(AF, dB_limit)
        plt.plot(x * (180 / np.pi), AF); plt.axis('tight'); plt.grid(True)        
        plt.title('Far-field Pattern')
        plt.xlabel(r'$\phi$')
        plt.ylabel('AF')
        plt.show()        
                     
    return m_index, zm

def FS_Chebyshev(N, R, x_min, x_max, x_num, plot_far=False, dB_limit= -40):
    """Function to evaluate Fourier series coefficients of Chebyshev far-field
       pattern"""
       
    c = np.cosh(np.arccosh(R) / (N))
    c = str(c)     
           
    if(N % 2 == 0):
        m_start = int(-N / 2)
        m_stop = int(N / 2)
        N = str(N)       
        fun_str_re = 'special.eval_chebyt(' + N + ',' + c + '*np.cos(x/2))'
        m_index, zm = FS(fun_str_re, m_start=m_start, m_stop=m_stop, err_lim=1e-5)        
    else:
        m_start = -N # make this (2*P+1) ... and take fourier for only half period
        m_stop = N
        N = str(N)
        fun_str_re = 'special.eval_chebyt(' + N + ',' + c + '*np.cos(x))'
        m_index, zm = FS(fun_str_re, m_start=m_start, m_stop=m_stop, err_lim=1e-5)
        
    if(plot_far):
        x, AF = IFS(zm, 2 * np.pi, m_start, m_stop, x_min, x_max, x_num)        
        AF = 20 * np.log10(abs(AF))
        AF = pl.cutoff(AF, dB_limit)
        plt.plot(x * (180 / np.pi), AF); plt.axis('tight'); plt.grid(True)        
        plt.title('Far-field Pattern')
        plt.xlabel(r'$\phi$')
        plt.ylabel('AF')
        plt.show()        
                     
    return m_index, zm

def FS_Zolotarev(N, R, x_min, x_max, x_num, plot_far=False, dB_limit= -40):
    """Function to evaluate Fourier series coefficients of Chebyshev far-field
       pattern"""    
           
    if(N % 2 == 0):
        print "Order needs to be an ODD number for null patterns"
    else:
        m_start = -N # make this (2*P+1) ... and take fourier for only half period
        m_stop = N
        m = zl.z_m_frm_R(N, R, a=0.1, b=0.9999999999999)
        m = str(m)        
        N = str(N)
        fun_str_re = 'zl.z_Zolotarev(' + N + ',' + 'np.sin(x)' + ',' + m +')'
        m_index, zm = FS(fun_str_re, m_start=m_start, m_stop=m_stop, err_lim=1e-5)
        
    if(plot_far):
        x, AF = IFS(zm, 2 * np.pi, m_start, m_stop, x_min, x_max, x_num)        
        AF = 20 * np.log10(abs(AF))
        AF = pl.cutoff(AF, dB_limit)
        plt.plot(x * (180 / np.pi), AF); plt.axis('tight'); plt.grid(True)        
        plt.title('Far-field Pattern')
        plt.xlabel(r'$\phi$')
        plt.ylabel('AF')
        plt.show()        
                     
    return m_index, zm

def dist_c_az(P, N, radius, R, mbar=5, alpha_x=0, starts_at_zero=True, 
              dist_type=None,plot_far=False, plot_modes=False, scan=False):
    r"""
    This function gives array excitation coefficients corresponding to various
    'circular' array distribution types (for 'azhimutal' patterns) such as 
    Chebyshev, Zolotarev, Taylor, Bayliss, etc.

    :param P:           order of the 'continuous' distribution (be careful! it 
                        has nothing to do with the number of elements 'N')
    :param N:           number of elements in the uniformly spaced circular ring 
                        array              
    :param radius:      radius of the circular ring in wavelengths   
    :param R:           side-lobe ratio in linear scale
    :param dist_type:   type of the distribution, e.g., 'Chebyshev', 'Zolotarev', etc
    :param bar:         transition index for dilation
    :param alpha:       Taylor's asymptotic tapering parameter
    
    :rtype:             A, a Numpy array of size (N, 1)
    """
        
    if(dist_type == 'Chebyshev'):
        m, zm = FS_Chebyshev(P, R, x_min=0, x_max=2 * np.pi, x_num=500, 
                             plot_far=plot_far, dB_limit= -40)
    elif(dist_type == 'Taylor'):
        m, zm = FS_Taylor(P, R, mbar, alpha_x, x_min=0, x_max=2*np.pi, 
                          x_num=500, plot_far=plot_far, dB_limit= -40)
    elif(dist_type == 'Zolotarev'):
        m, zm = FS_Zolotarev(P, R, x_min=0, x_max=2 * np.pi, x_num=500, 
                             plot_far=plot_far, dB_limit= -40)
    elif(dist_type == 'Bayliss'):
        m, zm = FS_Bayliss(P, R, mbar, alpha_x, x_min=0, x_max=2*np.pi, 
                          x_num=500, plot_far=plot_far, dB_limit= -40)        
        
    cm = zm / ((1j ** m) * special.jn(m, 2 * np.pi * radius))
        
    # Plotting the normalized absolute values of zm and cm
    if(plot_modes):
        plt.plot(m, abs(cm) / abs(cm).max(), 'r--', label=r"$\mathrm{a_m}$")
        plt.plot(m, abs(zm) / abs(zm).max(), label=r"$\mathrm{f_m}$")
        plt.axis('tight'); plt.grid(True)    
        plt.title(r'$\mathrm{Far\ \&\  Near-field\ \ Modes \ (f_m\ \& \ a_m)}$')
        plt.xlabel(r'$\mathrm{m\ (Mode\ Number)}$')
        plt.ylabel(r'$\mathrm{Mode\ Amplitude\ (abs)}$')
        plt.legend(bbox_to_anchor=(1, 1), loc=1, borderaxespad=1)
        plt.show()
        
    # Finally, discretizing the continuous distribution    
    if (starts_at_zero):
        beta = (np.linspace(1, N, num=N) - 1) * (2 * np.pi / N)
    else:
        beta = (np.linspace(1, N, num=N) - 0.5) * (2 * np.pi / N)
    beta = np.reshape(beta, (N, -1))    
    beta_tile = np.tile(beta, cm.size)        
    m_tile = np.tile(m.T, (N, 1))
    A = np.dot(np.exp(1j * m_tile * (beta_tile - scan)), cm)

    return A, cm

if __name__ == '__main__':
    

#==============================================================================
# Circular Chebyshev related script
#==============================================================================

    SLR = 25 # side-lobe ratio in dB
    R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
    P = 21
    N = 150
    radius = 2 / np.pi   
    
    A, cm = dist_c_az(P, N, radius, R, mbar=5, alpha_x=0, starts_at_zero=True, 
              dist_type='Chebyshev',plot_far=False, plot_modes=False, scan=False)
    array_ip = ip_format_c(N, radius, A, starts_at_zero=True, plot_type="3d", stem=False)

    phi, F = pl.pattern_t(array_ip, tht_scan=(0) * np.pi, phi_scan=(0) * np.pi, tht=0.5 * np.pi,
                 phi_min=0, phi_max=2 * np.pi, phi_num=500, scale="dB",
                 dB_limit= -200, factor="NF", plot_type='rect')
    
#==============================================================================
# Circular Taylor related script
#==============================================================================

#    SLR = 25 # side-lobe ratio in dB
#    R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
#    P = 9
#    N = 150
#    radius = 2 / np.pi
#    mbar = 3
#    alpha_x =0
#    
##    A = pl.dist(1, P + 1, R, dist_type_x='Taylor', mbar=mbar, alpha_x=alpha_x)
##    array_ip = pl.ip_format(1, 0, A, plot=True, stem=True, mayavi_app=False)
##    pl.pattern_u(array_ip, u_scan=0, u_min= 0, u_max=1, u_num=700, scale="dB",
##          dB_limit= -40, factor="AF", plot_type="rect", lattice=True)    
#    
##    print eval_Bayliss(P, R, mbar, alpha_x, x=0)
##    m_index, zm = FS_Bayliss(P, R, mbar, alpha_x, x_min=0, x_max=2*np.pi, x_num=500, plot_far=True, dB_limit= -40)
#    
#    A, cm = dist_c_az(P, N, radius, R, mbar, alpha_x, starts_at_zero=True, 
#              dist_type='Zolotarev',plot_far=True, plot_modes=False, scan=False)
#    array_ip = ip_format_c(N, radius, A, starts_at_zero=True, plot_type="3d", stem=False)
#
#    phi, F = pl.pattern_t(array_ip, tht_scan=(0) * np.pi, phi_scan=(0) * np.pi, tht=0.5 * np.pi,
#                 phi_min=0, phi_max=2 * np.pi, phi_num=500, scale="dB",
#                 dB_limit= -120, factor="NF", plot_type='polar')  

#==============================================================================
# Notes to myself
#============================================================================== 

# streamline Fourier coefficient evaluation ... i.e., FS_Chebyshev, etc 
    