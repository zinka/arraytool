#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

""" 
Module for analysis and design of planar array antennas.

Note:
Right now, this module contains only some important "basic routines".
In the near future, I will try to add more functions to this module.
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from enthought.mayavi import mlab

#==============================================================================
# Custom functions
#==============================================================================

def ip_format(a, b, A, gamma=np.pi / 2, plot=False, mayavi_app=False):
    """
    Function to generate the 'Arraytool' input format.
    
    a --> separation between elements along the x-axis in meters
    b --> separation between elements along the y-axis in meters
    A --> visual excitation matrix
    gamma --> lattice angle in radians
    plot --> if True, this function automatically produces a 3D plot. In 
             order to use this option you need MayaVi library.
    mayavi_app --> if True, 3D plot will be opened in the MayaVi main 
                   application itself
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
    array_ip = np.hstack((x, y, z, A))  # 'Arraytool' input format
    
    # plotting the array excitation, if plot option is True
    if (plot):
        if (mayavi_app): # this option opens the 3D plot in MayaVi Application
            mlab.options.backend = 'envisage'
        s1 = mlab.quiver3d(x, y, z, z, z, A) # stem3D representation
        ranges1 = [x.min(), x.max(), y.min(), y.max(), A.min(), A.max()]
        mlab.axes(xlabel="x", ylabel="y", zlabel="Excitation", ranges=ranges1)
        s1.scene.isometric_view()
        mlab.show()
        
    return array_ip

def pattern_u(array_ip, freq, u_scan=0, u_min= -1, u_max=1, u_num=50,
              plot=False, dB_limit= -40, gain=True, normalize=False,
              lattice=False):
    """
    Function to evaluate 2d array-factor(AF) or gain-factor(GF) of a 
    linear array in u-domain. By default, this function calculates gain-factor.
    
    array_ip --> input data in 'Arraytool' input format
    freq --> frequency of operation in Hzs
    u_scan --> beam scan position. For example, if you need to scan to 
               theta=30deg, u_scan = sin(pi/6). By default, 'u_scan' is 0.
    u_min, u_max --> limits of u-domain, by default they are -1 and +1 which 
                     correspond to the visible-space
    u_num --> number of sampling points between 'u_min' and 'u_max', including 
              boundaries. By default, is 50.
    plot --> By default, this is False. If 'plot="linear"', this function 
             produces a 2D plot in linear scale. Similarly, if 'plot="dB"',
             it produces 2D plot in dB scale. You need 'matplotlib' library
             to use this option.
    dB_limit --> when AF/GF/NF is 0, their dB value is -infinity. So, we will
                 restrict the minimum value to dB_limit. By default, 
                 'dB_limit' is -40.
    gain --> This function evaluates gain-factor(GF) by default. If you want,
             array-factor(AF), choose 'gain' as False. 
    normalize --> If 'normalize' is chosen as True, then both AF and GF 
                  become equal to normalized factor(NF) and possible 
                  maximum value of NF is 1 or 0dB. 
                  By default, 'normalize' is False.
    lattice --> by default, is False. If True, it will highlight both 
                visible-space and lattice period in u-domain.
    """
    
    x = array_ip[:, 0]
    y = array_ip[:, 1]
    z = array_ip[:, 2]
    A = array_ip[:, 3]
    k = (2 * np.pi * freq) / 3e8
    
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
        
        # Evaluating array-factor of the linear array
        AF = np.dot(np.exp(1j * k * U * X), A)
        
        # Evaluation of F=(AF or GF or NF) => depending upon the user's choice
        if(normalize):
            gain = False # because GF is not need anymore
            norm_fact = (abs(A)).sum()
            F = AF / norm_fact
            ff = "Gain or Array Factor "; f1 = "NF "; n1 = "Normalized "
        else:
            if(gain):            
                P_inc = ((abs(A)) ** 2).sum() # Converting the AF to GF
                GF = AF / np.sqrt(P_inc)
                F = GF; ff = "Gain-Factor "; f1 = "GF "; n1 = ""
            else:
                F = AF; ff = "Array-Factor "; f1 = "AF "; n1 = ""
                        
        # plotting the factor (AF/GF) in either linear or dB scale
        if(plot):
            if(plot == "linear"):
                plt.plot(u, F)
                plt.title(n1 + ff + "in linear scale")
            elif(plot == "dB"):
                F = 20 * np.log10(abs(F))
                # cutoff the "F" below some limit by using "masked_less"
                F = np.ma.masked_less(F, dB_limit)
                plt.plot(u, F) # use "F.data" for unmasked F
                plt.title(n1 + ff + "in dB scale")
            if(lattice):
                plt.axvspan(-1, +1, facecolor='y', alpha=0.2)
                a = x[2] - x[1]
                lim = -np.pi / (a * k)                
                plt.axvspan(-lim, +lim, facecolor='b', alpha=0.2)
            plt.axis('tight'); plt.grid(True)            
            plt.xlabel('u, where "u=sin(theta)" in visible-space')
            plt.ylabel(f1 + '(u)')
            plt.show()
                        
    return u, F

def AF_zeros(freq, a, M, R, type="DC"):
    """
    This function gives array-factor zeros corresponding to different 
    types of array distributions.
    """
    
    k = (2 * np.pi * freq) / 3e8
    m = np.ceil((M - 2) / 2)
    n = np.arange(1, 1 + m, 1)
    if(type == "DC"): # Dolph-Chebyshev zeros        
        c = np.cosh(np.arccosh(R) / (M - 1))
        K = (2 / a) * np.arccos((np.cos(np.pi * (2 * n - 1) / (2 * M - 2))) / c)
    elif(type == "RC"): # Riblet-Chebyshev zeros
        c1 = np.cosh(np.arccosh(R) / m)
        c = np.sqrt((1 + c1) / (2 + (c1 - 1) * np.cos(k * a / 2) ** 2))
        alph = c * np.cos(k * a / 2)
        xi = (1 / c) * np.sqrt(((1 + alph ** 2) / 2) + ((1 - alph ** 2) / 2) * 
                               np.cos(((2 * n - 1) * np.pi) / (2 * m)))
        K = (2 / a) * np.arccos(xi)
    elif(type=="MZs"): # McNamara-Zolotarev sum-pattern zeros
        K = "Yet to be done"
    elif(type=="MZd"): # McNamara-Zolotarev difference-pattern zeros
        K = "Yet to be done"
    K = np.reshape(K, (len(K), -1))        

    return K

def A_frm_zeros(K, a, M, symmetry=True):
        """
        My function description here.
        """
        
        sz = len(K)
        KK = np.tile(K, sz)
        if(symmetry):
            if(M % 2 == 0):
                tmp1 = np.arange(1, 1 + sz, 1) + 0.5
                tmp2 = 2 * np.cos(K * a / 2)
            else:
                tmp1 = np.arange(1, 1 + sz, 1)
                tmp2 = np.ones_like(K)
            tmp1 = np.reshape(tmp1, (-1, sz))
            TT = np.tile(tmp1, (sz, 1))
            CC = -np.linalg.inv(2 * np.cos(KK * TT * a))
            A = np.dot(CC, tmp2)
            A1 = np.flipud(A)
            
            if(M % 2 == 0):
                A_tot = np.vstack((A1, 1, 1, A))
            else:
                A_tot = np.vstack((A1, 1, A))
            
        return A_tot

#==============================================================================
# '__main__' function
#==============================================================================
    
if __name__ == '__main__':    

    #==========================================================================
    # frequency and array-arrangement
    #==========================================================================
    
    freq = 10e9 # frequency of operation in Hzs
    M = 7 # no. of elements along the x-axis
    N = 01 # no. of elements along the y-axis    
    a = 10e-3 # separation between elements along the x-axis in meters
    b = 10e-3 # separation between elements along the y-axis in meters
    gamma = np.pi / 2.5 # lattice angle in radians
    
    #==========================================================================
    # Array excitation
    #==========================================================================
    
    # Side-lobe information, if any exist!
    SLR = 30 # side-lobe ratio in dB
    R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
        
#    # Uniform
#    A = np.ones((N, M))
#    # Random
#    A = np.random.rand(N, M)

    # Dolph-Chebyshev type distributions    
    K = AF_zeros(freq, a, M, R, type="RC") # finding array-factor zeros    
    A = A_frm_zeros(K, a, M, symmetry=True).T # finding excitation coefficients
    
    #==========================================================================
    # Converting 'excitation & position' info into 'Arraytool' input format
    #==========================================================================
    
    array_ip = ip_format(a, b, A, gamma, plot=False, mayavi_app=False)
    
    #==========================================================================
    # Calling the 'pattern_u' function to evaluate and plot 2D AF/GF/NF
    #==========================================================================
    
    u_scan = 0
    u_min = -2
    u_max = +2
    u_num = 500     
    pattern_u(array_ip, freq, u_scan, u_min, u_max, u_num, plot="dB",
              dB_limit= -40, gain=True, normalize=True, lattice=True)
