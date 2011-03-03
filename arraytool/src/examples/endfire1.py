#! /usr/bin/env python
""" My test script. """

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
### for Palatino and other serif fonts use:
##rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'

from arraytool.planar import *

def pattern_u1(array_ip, u_scan=0, u_min= -1, u_max=1, u_num=50, scale="dB",
              dB_limit= -40, factor="GF", plot_type="rect", lattice=False):
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
        th = np.arcsin(u)
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
#            norm_fact = (abs(A)).sum()
            U2 = np.tile(1, len(x))
            norm_fact = np.dot(np.exp(1j * k * U2 * X), A)        
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
            if(plot_type == "rect"): # rectangular plot
                plt.plot(u, F_plt, 'k', label="$\mathrm{Uni}$ - $\mathrm{directional}$") # use "F.data" for unmasked F, if any exist                
                if(lattice): # highlighting visible-space and unit-lattice
                    plt.axvspan(-1, +1, facecolor='y', alpha=0.2)
                    lim = -np.pi / ((x[2] - x[1]) * k)                
                    plt.axvspan(-lim, +lim, facecolor='b', alpha=0.2)
                plt.axis('tight'); plt.grid(True)            
                plt.xlabel('u, where "u=sin(theta)" in visible-space')
                plt.ylabel(f1 + '(u)')                
            if(plot_type == "polar"): # polar plot
                if(scale == "linear"):
                    plt.polar(th, F_plt)
                    plt.polar(np.pi - th, F, '-b')
                if(scale == "dB"):
                    plt.polar(th, F_plt - dB_limit)
                    plt.polar(np.pi - th, F_plt - dB_limit, '-b')                
#            plt.title(n1 + ff + ss)
            plt.legend(bbox_to_anchor=(1.05, 1.05), loc=1, borderaxespad=1.)
#            plt.show()
                        
    return u, F

def pattern_u2(array_ip, u_scan=0, u_min= -1, u_max=1, u_num=50, scale="dB",
              dB_limit= -40, factor="GF", plot_type="rect", lattice=False):
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
        th = np.arcsin(u)
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
#            norm_fact = (abs(A)).sum()
            U2 = np.tile(1-u_scan, len(x))
            norm_fact = np.dot(np.exp(1j * k * U2 * X), A)        
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
            if(plot_type == "rect"): # rectangular plot
                p1 = plt.subplot(111)
                p1.plot(u, F_plt, '--k', label="$\mathrm{Bi}$ - $\mathrm{directional}$") # use "F.data" for unmasked F, if any exist                
                if(lattice): # highlighting visible-space and unit-lattice
                    plt.axvspan(-1, +1, facecolor='k', alpha=0.1)
                    lim = -np.pi / ((x[2] - x[1]) * k)                
                    plt.axvspan(-lim, +lim, facecolor='k', alpha=0.2)
                plt.axis('tight'); plt.grid(True)            
                plt.xlabel(r'$k_x/k_0$',fontsize=20)
                plt.ylabel(r'$\mathrm{Normalized}$ $AF$ ($k_x$), $\mathrm{in}$ $\mathrm{dB}$',fontsize=20)                
            if(plot_type == "polar"): # polar plot
                if(scale == "linear"):
                    plt.polar(th, F_plt)
                    plt.polar(np.pi - th, F, '-b')
                if(scale == "dB"):
                    plt.polar(th, F_plt - dB_limit)
                    plt.polar(np.pi - th, F_plt - dB_limit, '-b')                
#            plt.title(r'$\mathrm{Normalized}$ $\mathrm{Arrayfactor}$ $\mathrm{in}$ $\mathrm{dB}$ $\mathrm{scale}$',fontsize=20)
            #Text box
            props = dict(boxstyle='round', facecolor='w', alpha=0.5)
            textstr = '$M=%.0f$\n$d=%.2f\lambda$\n$\mathrm{SLR}=%.0f$$\mathrm{dB}$'%(M, a, 20)
            p1.text(0.14, 0.2, textstr, transform=p1.transAxes, fontsize=14,
                    verticalalignment='top', bbox=props)
            p1.legend(bbox_to_anchor=(0, 1), loc=2, borderaxespad=1)
            p1.axis(fontsize = 40)
                        
    return u, F

# actual values
freq = 10e9 # frequency of operation in Hzs
wav_len = 3e8 / freq # wavelength in meters
M = 6 # no. of elements along the x-axis
N = 5 # no. of elements along the y-axis       
a1 = 7.5e-3 # separation between elements along the x-axis in meters
b1 = 17e-3 # separation between elements along the y-axis in meters
gamma = np.pi / 2 # lattice angle in radians

# normalized values   
a = a1 / wav_len # 'a1' in-terms of lambda (wavelength)
b = b1 / wav_len # 'b1' in-terms of lambda (wavelength)

# Array excitation
SLR = 20 # side-lobe ratio in dB
R = 10 ** (SLR / 20) # converting SLR from dB scale to linear scale
#==============================================================================
# ARRAY 1
#==============================================================================
# Dolph-Chebyshev type distributions    
U0 = AF_zeros(a, M, R, type="DuhamelU") # finding array-factor zeros
print U0
A = A_frm_zeros(U0, a, M, symmetry=False).T # finding excitation coefficients
print A.T    

# Converting 'excitation & position' info into 'Arraytool' input format
array_ip = ip_format(a, b, A, gamma, plot=False, mayavi_app=False)

# Calling the 'pattern_u' function to evaluate and plot 2D AF/GF/NF
pattern_u1(array_ip, u_scan=0, u_min=-4, u_max=4, u_num=500, scale="dB",
      dB_limit= -40, factor="NF", plot_type="rect", lattice=False)
#==============================================================================
# ARRAY 2
#==============================================================================
# Dolph-Chebyshev type distributions    
U0 = AF_zeros(a, M, R, type="DuhamelB") # finding array-factor zeros
print U0
A = A_frm_zeros(U0, a, M, symmetry="odd").T # finding excitation coefficients
print A.T    

# Converting 'excitation & position' info into 'Arraytool' input format
array_ip = ip_format(a, b, A, gamma, plot=False, mayavi_app=False)

u_scan = 0
# Calling the 'pattern_u' function to evaluate and plot 2D AF/GF/NF
pattern_u2(array_ip, u_scan=u_scan, u_min=-4, u_max=4, u_num=1000, scale="dB",
      dB_limit= -40, factor="NF", plot_type="rect", lattice=True)

plt.show()