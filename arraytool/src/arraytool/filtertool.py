#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2014 Srinivasa Rao Zinka
# License: New BSD License.

from __future__ import division
import numpy as np
import sympy as sp
import scipy.signal as signal
import matplotlib.pyplot as plt
import warnings
import Tkinter as ti
import tkFileDialog as tkdlg

# adjusting "matplotlib" label fonts
from matplotlib import rc
rc('text', usetex=True)

def ft_import(dtype='complex'):
    r"""
    A simple function to import a CSV text file as a Numpy ndarray
    
    :param dtype: Data-type of the resulting array; default:complex. For further
                  information, see numpy.loadtxt.
                  
    :rtype:       Returns a Numpy ndarray
    """
    master = ti.Tk(); master.withdraw() #hiding tkinter window 
    file_path = tkdlg.askopenfilename(title="Open file", filetypes=[("txt file",
                ".csv"), ("All files", ".*")]); master.quit()
    ip = np.loadtxt(file_path, delimiter=',', dtype=dtype)
    return ip

def cutoff(F, dB_limit= -40):
    r"""
    When magnitude of S11 or S21 is 0 in linear scale, their dB value is '-infinity'.
    So, this function will be used to cut-off all the value below some 'dB_limit'.
    
    :param F:           F is a Numpy array
    :param dB_limit:    cut-off level in dB, default value is -40
                  
    :rtype:             a Numpy array, same size as the input array F
    """
    msk1 = F < dB_limit
    fill = msk1 * dB_limit
    msk2 = F >= dB_limit
    F = F * (msk2) + fill
    return F

def I_to_i(A):
    r"""
    Function to convert the symbol 'I' of SymPy to 'i'.
    
    :param A: 
                  
    :rtype:    
    """
    B = np.zeros((len(A), 1), dtype=complex)
    for i in range(len(A)):
        B[i] = complex(A[i])
    return B

def s_to_w(poly_ip, coef_norm=False):
    r"""
    Arraytool mostly uses polynomials defined in 's' domain, i.e., the Laplace domain.
    However, sometimes we may need polynomials in 'w' (omega) domain, i.e., lowpass
    prototype frequency domain. So, this function can be used to convert a given
    polynomial from 's' to 'w' domain (s=jw).
    
    :param poly_ip:    Input polynomial
    :param coef_norm:  If `true`, the output polynomial will be normalized such
                       that the coefficient of the highest degree term becomes 1.
                  
    :rtype:            Returns a polynomial of the same order, however defined
                       in `s=jw` domain
    """
    if(len(poly_ip) == 1):
        poly_op = np.array([[1]])
    else:        
        s = sp.Symbol('s'); w = sp.Symbol('w')
        poly_ip = sp.Poly(poly_ip.ravel().tolist(), s)
        poly_op = sp.simplify(poly_ip.subs(s, 1j * w)) # substitution
        poly_op = sp.Poly(poly_op, w).all_coeffs()
        poly_op = I_to_i(poly_op)
        poly_op = np.reshape(poly_op, (len(poly_op), -1))
        if(coef_norm): poly_op = poly_op / poly_op[0]
    return poly_op

def w_to_s(poly_ip, coef_norm=False):
    r"""
    Arraytool mostly uses polynomials defined in 's' domain, i.e., Laplace domain.
    So, this function can be used to convert a given polynomial from 'w' to 's' 
    domain (w=-js).
    
    :param poly_ip:     Input polynomial
    :param coef_norm:   If `true`, the output polynomial will be normalized such
                        that the coefficient of the highest degree term becomes 1.
                  
    :rtype:             Returns a polynomial of the same order, however defined
                        in `w=-js` domain
    """
    if(len(poly_ip) == 1):
        poly_op = np.array([[1]])
    else:
        s = sp.Symbol('s'); w = sp.Symbol('w')
        poly_ip = sp.Poly(poly_ip.ravel().tolist(), w)
        poly_op = sp.simplify(poly_ip.subs(w, -1j * s)) # substitution
        poly_op = sp.Poly(poly_op, s).all_coeffs()
        poly_op = I_to_i(poly_op)
        poly_op = np.reshape(poly_op, (len(poly_op), -1))
        if(coef_norm): poly_op = poly_op / poly_op[0]
    return poly_op

def plot_rational(Num, Den, x_min= -1, x_max=1, x_num=100):
    r"""
    Simple function to plot a given rational function (i.e., a ratio of two 
    polynomial functions) of x.
       
    :param Num:    Numerator polynomial as a function of x
    :param Den:    Denominator polynomial as a function of x
    :param x_min:  lower limit of x (for plotting)  
    :param x_max:  upper limit of x (for plotting) 
    :param x_num:  number of points between 'x_min' and 'x_max' including x_min and x_max  
                  
    :rtype:        [x,y] ... one can simply use these parameter as inputs to plot(x,y)
    """
    x = np.linspace(x_min, x_max, x_num)
    y_Num = np.polyval(Num, x); y_Den = np.polyval(Den, x); y = y_Num / y_Den    
    plt.plot(x, y); plt.axis('tight'); plt.grid(True)
    plt.xlabel('x'); plt.ylabel('y'); plt.show()
    return x, y

def Gramm_Schmidt(A):
    r"""
    Function to implement the modified Gram-Schmidt orthonormalization.
    
    :param A:  A is the "partial" coupling matrix, where we know only first and last rows.
    
    :rtype:    Returns the final coupling matrix using Gramm-Schmidt method  
    """
    A = A.T; m, n = A.shape
    q = np.zeros((m, n), dtype='complex'); r = np.zeros((n, n), dtype='complex')
    for k in range(n):
        r[k, k] = np.linalg.norm(A[:, k])
        q[:, k] = A[:, k] / r[k, k]
        r[k, k + 1:n] = np.dot(q[:, k], A[:, k + 1:n])
        A[:, k + 1:n] = A[:, k + 1:n] - np.outer(q[:, k], r[k, k + 1:n])
    return q.T

def Chebyshev_gen(N, poles):
    r"""
    Function to evaluate the numerator and denominator polynomials of the 
    generalized Chebyshev (i.e., with/without poles) filtering function.
    
    :param N:       Order of the generalized Chebyshev polynomial
    :param poles:   np.array([pole1, pole2, ...]) ... where pole1, pole2, etc,
                    are poles of the generalized Chebyshev polynomial located at 
                    finite frequencies, i.e., not at w = infinity.
                  
    :rtype:         Returns numerator and denominator as polynomials
    """
    if(len(poles) == 0):
        Num = sp.polys.orthopolys.chebyshevt_poly(N)
        Num = sp.Poly(Num).all_coeffs()
        Den = np.array([[1]])
    else:        
        # evaluation of the polynomial P (i.e., the denominator)
        P = np.poly(poles); Den = np.reshape(P, (len(P), -1)) # denominator
        # placing all other poles at infinity
        poles = np.hstack((poles, np.inf * np.ones(N - poles.size)))
        # R. J. Cameron's recursive algorithm
        x = sp.Symbol('x')
        tmp1 = x - 1 / poles[0]
        tmp2 = sp.sqrt(x ** 2 - 1) * np.sqrt(1 - 1 / poles[0] ** 2)
        for j in range(1, len(poles)):
            if(np.isinf(poles[j])):
                U = x * tmp1 + sp.sqrt(x ** 2 - 1) * tmp2
                V = x * tmp2 + sp.sqrt(x ** 2 - 1) * tmp1
                tmp1 = sp.simplify(U); tmp2 = sp.simplify(V)
            else:
                U = x * tmp1 - tmp1 / poles[j] + sp.sqrt(x ** 2 - 1) * np.sqrt(1 - 1 / poles[j] ** 2) * tmp2
                V = x * tmp2 - tmp2 / poles[j] + sp.sqrt(x ** 2 - 1) * np.sqrt(1 - 1 / poles[j] ** 2) * tmp1
                tmp1 = sp.simplify(U); tmp2 = sp.simplify(V)
        U = tmp1; Num = sp.Poly(U, x).all_coeffs() # numerator
        Num = I_to_i(Num)
    norm = np.polyval(Num, 1) / np.polyval(Den, 1) 
    Num = np.reshape(Num, (len(Num), -1)) / norm # so that abs(polynomial) value becomes 1 at +1
    return Num, Den

def poly_E(eps, eps_R, F, P):
    r"""
    Function to obtain the polynomial E and its roots in the s-domain.
    
    :param eps:     Constant term associated with S21
    :param eps_R:   Constant term associated with S11
    :param F:       Polynomial F, i.e., numerator of S11 (in s-domain)
    :param P:       Polynomial P, i.e., numerator of S21 (in s-domain)
    
    For further explanation, see "filter theory notes" by the author.
                  
    :rtype:         [polynomial_E, roots_E] ... where polynomial E is the 
                    denominator of S11 and S21
    """
    N = len(F); nfz = len(P)
    if(((N+nfz)%2==0) and (abs(eps.real)>0)):
        warnings.warn("'eps' value should be pure imaginary when (N+nfz) is an even number")
    elif(((N+nfz+1)%2==0) and (abs(eps.imag)>0)):
        warnings.warn("'eps' value should be pure real when (N+nfz) is an odd number")
    s = sp.Symbol('s')
    poly_P = sp.Poly(P.ravel().tolist(), s)
    poly_F = sp.Poly(F.ravel().tolist(), s)
    poly_E = eps_R * poly_P + eps * poly_F
    roots_E = I_to_i(sp.nroots(poly_E))
    roots_E = np.reshape(roots_E, (len(roots_E), -1))
    roots_E = -abs(roots_E.real) + 1j * roots_E.imag # all roots to the LHS
    # create the polynomial from the obtained LHS roots
    poly_E = np.poly(roots_E.ravel().tolist())
    poly_E = np.reshape(poly_E, (len(poly_E), -1))
    return poly_E, roots_E

def plot_mag(eps, eps_R, F, P, E, w_min= -2, w_max=2, w_num=500, dB=True,
             dB_limit= -40, plot=True):
    r"""
    Function to plot magnitudes of S11 and S21 in either linear or dB scale.
    
    :param eps:     Constant term associated with S21
    :param eps_R:   Constant term associated with S11
    :param F:       Polynomial F, i.e., numerator of S11 (in s-domain)
    :param P:       Polynomial P, i.e., numerator of S21 (in s-domain)
    :param E:       polynomial E is the denominator of S11 and S21 (in s-domain)
    :param w_min:   lower limit of w (for plotting)  
    :param w_max:   upper limit of w (for plotting)  
    :param w_num:   number of points between 'x_min' and 'x_max' including x_min and x_max  
    :param dB:      If true plotting will be done in dB scale
    :param dB_limit: cut-off level in dB, default value is -40
    :param plot:     If True, plot will be shown
                  
    :rtype:          [w, S11, S21] ... frequency and the corresponding S11 and S21
                     all in linear scale
    """
    N = len(F); nfz = len(P)
    if(((N+nfz)%2==0) and (abs(eps.real)>0)):
        warnings.warn("'eps' value should be pure imaginary when (N+nfz) is an even number")
    elif(((N+nfz+1)%2==0) and (abs(eps.imag)>0)):
        warnings.warn("'eps' value should be pure real when (N+nfz) is an odd number")    
    w = np.linspace(w_min, w_max, w_num)
    F = s_to_w(F); P = s_to_w(P); E = s_to_w(E)
    F_val = np.polyval(F, w); P_val = np.polyval(P, w); E_val = np.polyval(E, w)
    S11 = (1 / eps_R) * (F_val / E_val); S21 = (1 / eps) * (P_val / E_val)
    if(plot):
        if(dB):
            S11_plt = 20 * np.log10(abs(S11)); S21_plt = 20 * np.log10(abs(S21))
            S11_plt = cutoff(S11_plt, dB_limit); S21_plt = cutoff(S21_plt, dB_limit)
            y_labl = r'$\ \mathrm{(dB)}$'
        else:
            S11_plt = abs(S11); S21_plt = abs(S21)
            y_labl = r'$\ \mathrm{(linear)}$'
        plt.plot(w, S21_plt, 'b-', label=r"$S_{21}$")
        plt.plot(w, S11_plt, 'r-', label=r"$S_{11}$")
        plt.axis('tight'); plt.grid(True); plt.legend()
        plt.xlabel(r'$\Omega\ \mathrm{(rad/s)}$', fontsize=14)
        plt.ylabel(r'$\mathrm{Magnitude}$' + y_labl, fontsize=14)    
        plt.show()
    
    return w, S11, S21

def plot_delay(roots_E, w_min= -2, w_max=2, w_num=500, plot=True):
    r"""
    Function to plot the group delay.
    
    :param roots_E:   Roots of the polynomial E (in s-domain). These roots can be either entered 
                      manually or can be obtained using the function `poly_E`.
    :param w_min:     lower limit of w (for plotting)  
    :param w_max:     upper limit of w (for plotting)
    :param w_num:     number of points between 'x_min' and 'x_max' including x_min and x_max
    :param plot:      If True, plot will be shown
                  
    :rtype:           Return [w,tau] ... where w is frequency and `tau` is group delay.
    """
    w = np.linspace(w_min, w_max, w_num)
    tau = np.zeros_like(w)
    for i in range(len(w)):
        tmp = 0
        for j in range(len(roots_E)):
            tmp += -roots_E[j].real / (roots_E[j].real ** 2 + (w[i] - roots_E[j].imag) ** 2)
        tau[i] = tmp
    if(plot):
        plt.plot(w, tau)
        plt.axis('tight'); plt.grid(True)
        plt.xlabel(r'$\Omega\ \mathrm{(rad/s)}$', fontsize=14)
        plt.ylabel(r'$\mathrm{Group}\ \mathrm{delay} \mathrm{\ (s)}$', fontsize=14)    
        plt.show()        
    return w, tau

def coupling_N(F, P, E, eps, eps_R):
    r"""
    Function to evaluate the (N,N) coupling matrix.
    
    :param F:       Polynomial F, i.e., numerator of S11 (in s-domain)
    :param P:       Polynomial P, i.e., numerator of S21 (in s-domain)
    :param E:       polynomial E is the denominator of S11 and S21 (in s-domain)
    :param eps:     Constant term associated with S21
    :param eps_R:   Constant term associated with S11
                  
    :rtype:         [M, RS_L1, RL_LN] ... M is the NbyN coupling matrix and RS_L1 and RL_LN
                    are the ratios of RS/L1 and RL/LN.
    """
    N = len(F); nfz = len(P)
    if(((N+nfz)%2==0) and (abs(eps.real)>0)):
        warnings.warn("'eps' value should be pure imaginary when (N+nfz) is an even number")
    elif(((N+nfz+1)%2==0) and (abs(eps.imag)>0)):
        warnings.warn("'eps' value should be pure real when (N+nfz) is an odd number")    
    F = s_to_w(F); P = s_to_w(P); E = s_to_w(E)
    nfz = len(P) - 1
    const_mult = np.conjugate(eps) / eps * (-1) ** nfz
    EF_plus = E + F / eps_R
    EF_plus_conj = const_mult * EF_plus.conj()
    EF_minus = E - F / eps_R
    EF_minus_conj = const_mult * EF_minus.conj()
    y11_Num = 1j * (EF_minus + EF_minus_conj)
    y21_Num = 2j * P / eps
    y_Den = EF_plus - EF_plus_conj 
    # The function "signal.residue" takes only 1D arrays!
    resid11, poles11, const11 = signal.residue(y11_Num[:, 0], y_Den[:, 0])
    resid21, poles21, const21 = signal.residue(y21_Num[:, 0], y_Den[:, 0])
    # Gramm_Schmidt orthonormalization
    T1k = np.sqrt(resid11); lambdk = -poles11; TNk = resid21 / T1k
    RS_L1 = sum(T1k ** 2); RL_LN = sum(TNk ** 2)    
    T = np.eye(len(T1k), len(T1k), dtype='complex')
    T[0, :] = T1k; T[1, :] = TNk
    np.set_printoptions(precision=6, suppress=True)
    T = Gramm_Schmidt(T) # "normalizing of T1k, TNk" is done in this step
    # swapping the second and last rows after normalization is finished
    temp = np.copy(T[1, :]); T[1, :] = T[-1, :]; T[-1, :] = temp
    Lamb = np.diag(lambdk) # diagonal eigenvalue matrix
    M = np.dot(T, np.dot(Lamb, T.T)) # (N,N) coupling matrix
    return M, RS_L1, RL_LN

def MN_to_Sparam(M, Rs, Rl, w_min= -2, w_max=2, w_num=500, dB=True,
                 dB_limit= -40, plot=True):
    r"""
    Function to plot S parameters from a given (N,N) coupling matrix.
    
    :param M:      NbyN coupling matrix
    :param Rs:     Source resistance
    :param Rl:     Load resistance
    :param w_min:  lower limit of w (for plotting) 
    :param w_max:  upper limit of w (for plotting) 
    :param w_num:  number of points between 'x_min' and 'x_max' including x_min and x_max
    :param dB:     If true plotting will be done in dB scale
    :param dB_limit:  cut-off level in dB, default value is -40
    :param plot:   If True, plot will be drawn ... but will be showed only if show = True
                  
    :rtype:        [w, S11, S21] ... frequency and the corresponding S11 and S21
                   all in linear scale   
    """
    w = np.linspace(w_min, w_max, w_num)
    R = np.zeros_like(M); R[0, 0] = Rs; R[-1, -1] = Rl
    MR = M - 1j * R ; I = np.eye(M.shape[0], M.shape[1])
    # Calculating S parameters
    S11 = np.zeros((len(w), 1), dtype=complex)
    S21 = np.zeros((len(w), 1), dtype=complex) # 'dtype' is important
    for i in range(len(w)):
        A = MR + w[i] * I
        A_inv = np.linalg.inv(A)
        S11[i] = 1 + 2j * Rs * A_inv[0, 0]
        S21[i] = -2j * np.sqrt(Rs * Rl) * A_inv[-1, 0]
    if(plot): # Plotting     
        # Converting the S parameters into either linear or dB scale
        if(dB):
            S11_plt = 20 * np.log10(abs(S11)); S21_plt = 20 * np.log10(abs(S21))
            S11_plt = cutoff(S11_plt, dB_limit); S21_plt = cutoff(S21_plt, dB_limit)
            y_labl = r'$\ \mathrm{(dB)}$'
        else:
            S11_plt = abs(S11); S21_plt = abs(S21)
            y_labl = r'$\ \mathrm{(linear)}$'
        plt.plot(w, S21_plt, 'b-', label=r"$S_{21}$")
        plt.plot(w, S11_plt, 'r-', label=r"$S_{11}$")
        plt.axis('tight'); plt.grid(True); plt.legend()
        plt.xlabel(r'$\Omega\ \mathrm{(rad/s)}$', fontsize=14)
        plt.ylabel(r'$\mathrm{Magnitude}$' + y_labl, fontsize=14)
        plt.show()        
    return w, S11, S21

def coupling_N2(F, P, E, eps, eps_R):
    r"""
    Function to evaluate the (N+2,N+2) coupling matrix.
    
    :param F:       Polynomial F, i.e., numerator of S11 (in s-domain)
    :param P:       Polynomial P, i.e., numerator of S21 (in s-domain)
    :param E:       polynomial E is the denominator of S11 and S21 (in s-domain)
    :param eps:     Constant term associated with S21
    :param eps_R:   Constant term associated with S11
                  
    :rtype:         Returns (N+2,N+2) coupling matrix 
    """
    F = s_to_w(F); P = s_to_w(P); E = s_to_w(E)
    nfz = len(P) - 1; N = len(E) - 1
    const_mult = np.conjugate(eps) / eps * (-1) ** nfz
    EF_plus = E + F / eps_R
    EF_plus_conj = const_mult * EF_plus.conj()
    EF_minus = E - F / eps_R
    EF_minus_conj = const_mult * EF_minus.conj()
    y11_Num = 1j * (EF_minus + EF_minus_conj)
    y21_Num = -2j * P / eps
    y_Den = EF_plus - EF_plus_conj
    # The function "signal.residue" takes only 1D arrays!
    resid11, poles11, const11 = signal.residue(y11_Num[:, 0], y_Den[:, 0])
    resid21, poles21, const21 = signal.residue(y21_Num[:, 0], y_Den[:, 0])
    MSk = np.sqrt(resid11); lambdk = -poles11; MLk = resid21 / MSk
    JSL = -const21
    M = np.zeros((N + 2, N + 2), dtype=complex)
    M[0, 1:N + 1] = np.reshape(MSk, (-1, N))
    M[-1, 1:N + 1] = np.reshape(MLk, (-1, N))
    M[:, 0] = M[0, :]; M[:, -1] = M[-1, :]
    M[0, -1] = JSL[0]; M[-1, 0] = JSL[0]
    diag1 = np.diag(lambdk); M[1:N + 1, 1:N + 1] = diag1
    return M

def MN2_to_Sparam(M, Rs=1, Rl=1, w_min= -2, w_max=2, w_num=500, dB=True,
                 dB_limit= -40, plot=True):
    r"""
    Function to plot S parameters from a given (N+2,N+2) coupling matrix.
    
    :param M:      NbyN coupling matrix
    :param Rs:     Source resistance
    :param Rl:     Load resistance
    :param w_min:  lower limit of w (for plotting) 
    :param w_max:  upper limit of w (for plotting) 
    :param w_num:  number of points between 'x_min' and 'x_max' including x_min and x_max
    :param dB:     If true plotting will be done in dB scale
    :param dB_limit:  cut-off level in dB, default value is -40
    :param plot:   If True, plot will be drawn ... but will be showed only if show = True
                  
    :rtype:    [w, S11, S21] ... frequency and the corresponding S11 and S21
               all in linear scale
    """
    w = np.linspace(w_min, w_max, w_num)
    R = np.zeros_like(M); R[0, 0] = Rs; R[-1, -1] = Rl
    MR = M - 1j * R ; I = np.eye(M.shape[0], M.shape[1])
    I[0, 0] = 0; I[-1, -1] = 0
    # Calculating S parameters
    S11 = np.zeros((len(w), 1), dtype=complex)
    S21 = np.zeros((len(w), 1), dtype=complex) # 'dtype' is important
    for i in range(len(w)):
        A = MR + w[i] * I
        A_inv = np.linalg.inv(A)
        S11[i] = 1 + 2j * Rs * A_inv[0, 0]
        S21[i] = -2j * np.sqrt(Rs * Rl) * A_inv[-1, 0]
    if(plot): # Plotting     
        # Converting the S parameters into either linear or dB scale
        if(dB):
            S11_plt = 20 * np.log10(abs(S11)); S21_plt = 20 * np.log10(abs(S21))
            S11_plt = cutoff(S11_plt, dB_limit); S21_plt = cutoff(S21_plt, dB_limit)
            y_labl = r'$\ \mathrm{(dB)}$'
        else:
            S11_plt = abs(S11); S21_plt = abs(S21)
            y_labl = r'$\ \mathrm{(linear)}$'
        plt.plot(w, S21_plt, 'b-', label=r"$S_{21}$")
        plt.plot(w, S11_plt, 'r-', label=r"$S_{11}$")
        plt.axis('tight'); plt.grid(True); plt.legend()
        plt.xlabel(r'$\Omega\ \mathrm{(rad/s)}$', fontsize=14)
        plt.ylabel(r'$\mathrm{Magnitude}$' + y_labl, fontsize=14)
        plt.show()        
    return w, S11, S21

def rotation_R(N,i,j,tht):
    r"""
    Function to calculate rotation matrix.
    
    :param N:   N is the order of the coupling matrix
    :param i:   One element of the pivot (i,j)
    :param j:   Another element of the pivot (i,j)
    :param tht: Rotation angle given in radians
                  
    :rtype:    Returns the corresponding rotation matrix of size (N,N)
    """
    R = np.eye(N)
    R[i-1,i-1]=np.cos(tht); R[j-1,j-1]=np.cos(tht)
    R[i-1,j-1]=-np.sin(tht); R[j-1,i-1]=np.sin(tht)
    return R

def  annihilation(p,q,i,j,N,M):
    r"""
    Function to annihilate a given element without changing eigen values.
    
    :param p:  Row index of the element to be annihilated
    :param q:  Column index of the element to be annihilated
    :param i:  One element of the pivot (i,j)
    :param j:  Another element of the pivot (i,j)
    :param N:  Order of the filter
    :param M:  M is the coupling matrix, either (N,N) or (N+2,N+2) matrix
                  
    :rtype:    Returns the rotation angle required to annihilate (p,q)th element
    """    
    if(p==q):
        if(p==i):
            print p==i
            tht = np.arctan((-M[i-1,j-1]+np.sqrt(M[i-1,j-1]**2-M[i-1,i-1]*M[j-1,j-1])/M[j-1,j-1]))
        else:
            print'p!=i'
            tht =np.arctan((M[i-1,j-1]+np.sqrt((M[i-1,j-1]**2-M[i-1,i-1]*M[j-1,j-1])))/M[i-1,i-1])            
    else:
        print p!=i
        if(p==i and q==j):
            tht=0.5*np.arctan(2*M[i-1,j-1]/(M[j-1,j-1]-M[i-1,i-1]))            
        elif (p==j and q==i):
            k=p            
            tht= 0.5*np.arctan((2*M[k-1,i-1])/(M[k-1,k-1]-M[j-1,j-1]))
        elif(p==i and q!=i and p!=j):
            k=p
            tht=np.arctan(M[k-1,q-1]/M[k,q-1])           
        elif(p==j and q!=i and q!=j):
            k = q            
            tht = -np.arctan(M[j-1,k-1]/M[i-1,k-1])
        elif(q==i and p!=i and q!=j):
            k = p            
            tht  = np.arctan(M[i-1,k-1]/M[j-1,k-1])            
        elif((q==j) and (p!=i) and (p!=j)):
            k = p
            tht  = -np.arctan(M[k-1,j-1]/M[k-1,i-1])           
        return tht
    
def K(plot=True):
    r"""
    Function to import (gap, f_even, f_odd) information as a CSV file and to
    evaluate coupling coefficients and plot them (vs gap), if necessary.
    
    :param plot:   If True, plot will be drawn
                  
    :rtype:        [gap,K] ... Returns gap and coupling coefficients
    """
    ip = ft_import(dtype='float')
    result = (ip[:,2]**2-ip[:,1]**2)/(ip[:,2]**2+ip[:,1]**2)
    if(plot):
        plt.plot(ip[:,0], result)
        plt.axis('tight'); plt.grid(True)
        plt.xlabel('gap')
        plt.ylabel('K')
        plt.show()
    return ip[:,0], result

def Lee_roots_map(roots, x1):
    r"""
    My function description here.
    
    :param roots: 
    :param x1:
                  
    :rtype:    
    """
    
    c1 = (1-x1); c2 = x1/(1-x1)
    roots1_p = 0.5*(c1*roots+np.sqrt(c1**2*roots**2-4*c1*c2))
    roots1_m = 0.5*(c1*roots-np.sqrt(c1**2*roots**2-4*c1*c2))
    roots1 = []
    
    for i in range(len(roots1_p)):
        if(roots1_p[i].imag >= 0):
            roots1.append(roots1_p[i])
            
    for i in range(len(roots1_m)):
        if(roots1_m[i].imag >= 0):
            roots1.append(roots1_m[i])
    roots1 = np.array(roots1)
    roots1 = np.concatenate((roots1, np.conjugate(roots1)))
    
    return roots1

def dual_band(F, P, E, eps, eps_R=1, x1=0.5, map_type=1):
    r"""
    Function to give modified F, P, and E polynomials after lowpass 
    prototype dual band transformation.
    
    :param F:       Polynomial F, i.e., numerator of S11 (in s-domain)
    :param P:       Polynomial P, i.e., numerator of S21 (in s-domain)
    :param E:       polynomial E is the denominator of S11 and S21 (in s-domain)
    :param eps:     Constant term associated with S21
    :param eps_R:   Constant term associated with S11
    :param x1:
    :param map_type:
                  
    :rtype:    
    """
    
    if((map_type == 1) or (map_type == 2)):
                
        s = sp.Symbol('s')        
        if (map_type == 1):
            a = -2j / (1 - x1 ** 2); b = -1j * (1 + x1 ** 2) / (1 - x1 ** 2)
        elif (map_type == 2):
            a = 2j / (1 - x1 ** 2); b = 1j * (1 + x1 ** 2) / (1 - x1 ** 2)        
        s1 = a * s ** 2 + b
                
        F = sp.Poly(F.ravel().tolist(), s)    
        F1 = sp.simplify(F.subs(s, s1))    
        F1 = sp.Poly(F1, s).all_coeffs()    
        F1 = I_to_i(F1)
        
        E = sp.Poly(E.ravel().tolist(), s)    
        E1 = sp.simplify(E.subs(s, s1))    
        E1 = sp.Poly(E1, s).all_coeffs()    
        E1 = I_to_i(E1)
        
        P = sp.Poly(P.ravel().tolist(), s)    
        P1 = sp.simplify(P.subs(s, s1))    
        P1 = sp.Poly(P1, s).all_coeffs()    
        P1 = I_to_i(P1)
        
    elif (map_type == 3):
        
        F_roots = np.roots(F.ravel().tolist())
        P_roots = np.roots(P.ravel().tolist())
        F1_roots = Lee_roots_map(F_roots, x1)
        P1_roots = Lee_roots_map(P_roots, x1)
        P1_roots = np.concatenate((P1_roots, np.array([0,0])))
        F1 = np.poly(F1_roots)
        P1 = np.poly(P1_roots)
        F1 = np.reshape(F1, (len(F1),-1))
        P1 = np.reshape(P1, (len(P1),-1))
        E1 = poly_E(eps, eps_R, F1, P1)[0]

    return F1, P1, E1

if __name__ == '__main__':
        
#==============================================================================
# 7th order example (P. 300, Sec. 8.3.1, R. J. Cameron et al.)
#==============================================================================    
    
    N = 7
    poles = np.array([1.2576, -0.1546 - 0.9218j, -0.1546 + 0.9218j])
#     poles = np.array([-3,3])
    eps = 6.0251j; eps_R = 1
#    poles = np.array([])
    F, P = Chebyshev_gen(N, poles)
    plot_rational(F, P, x_min= -1.1, x_max=1.1, x_num=1000)
 
    F = w_to_s(F, coef_norm=True)
    P = w_to_s(P, coef_norm=True) 
    
    print 'F:', '\n', F; print 'P:', '\n', P
    E, roots_E = poly_E(eps, eps_R, F, P)
    print 'E:', '\n', E
     
#    F1, P1, E1 = dual_band(F, P, E, 0.4)
     
    plot_mag(eps, eps_R, F, P, E, w_min= -2, w_max=2, w_num=500, dB=True,
             dB_limit= -40, plot=True)
    plot_delay(roots_E)
   
    # From now onwards, unlike the Cameron's example, this filter is doubly terminated
    M, Rs, Rl = coupling_N(F, P, E, eps, eps_R)
    print 'M:', '\n', M.real
    print 'Rs:', Rs
    print 'Rl:', Rl
     
    MN_to_Sparam(M, Rs, Rl, w_min= -3, w_max=3, w_num=500, dB=True, dB_limit= -40)
    
#==============================================================================
# 4th order example (P. 228, Sec. 6.3.2, R. J. Cameron et al.)
#==============================================================================

#    N = 4
#    poles = np.array([1.3217, 1.8082])
#    eps = 1.1548j; eps_R = 1
##    poles = np.array([])
#    F, P = Chebyshev_gen(N, poles)
##    plot_rational(F, P, x_min= -1, x_max=1, x_num=1000)
#    F = w_to_s(F, coef_norm=True)
#    P = w_to_s(P, coef_norm=True)    
#    print 'F:', '\n', F; print 'P:', '\n', P
#    print 'E:', '\n', poly_E(eps, eps_R, F, P)[0]
#    E = poly_E(eps, eps_R, F, P)[0]
#    plot_mag(eps, eps_R, F, P, E)

#==============================================================================
# 4th order example (P. 312, Sec. 8.4.2, R. J. Cameron et al.)
#==============================================================================

#    N = 4
#    poles = np.array([-3.7431, -1.8051, 1.5699, 6.1910])
#    eps = 33.140652j; eps_R = +1.000456
##    poles = np.array([])
#    F, P = Chebyshev_gen(N, poles)
##    plot_rational(F, P, x_min= -1, x_max=1, x_num=1000)
#    F = w_to_s(F, coef_norm=True)
#    P = w_to_s(P, coef_norm=True)    
#    print 'F:', '\n', F; print 'P:', '\n', P
#    [E, roots_E] = poly_E(eps, eps_R, F, P)
#    print 'E:', '\n', E
#    
#    np.set_printoptions(precision=4, suppress=True)    
#    M = coupling_N2(F, P, E, eps, eps_R)
#    print 'M:', '\n', M.real
#    MN2_to_Sparam(M, Rs=1, Rl=1, w_min= -8, w_max=8, w_num=500, dB=True, dB_limit= -50)
    
#    N = 4
#    poles = np.array([-6, -2, 2, 6])
#    eps = 33.140652j; eps_R = +1.000456
##    poles = np.array([])
#    F, P = Chebyshev_gen(N, poles)
##    plot_rational(F, P, x_min= -1, x_max=1, x_num=1000)
#    F = w_to_s(F, coef_norm=True)
#    P = w_to_s(P, coef_norm=True)    
#    print 'F:', '\n', F; print 'P:', '\n', P
#    [E, roots_E] = poly_E(eps, eps_R, F, P)
#    print 'E:', '\n', E
#    print E.ravel().tolist()
#    
#    np.set_printoptions(precision=4, suppress=True)    
#    M = coupling_N2(F, P, E, eps, eps_R)
#    print 'M:', '\n', M.real
#    MN2_to_Sparam(M, Rs=1, Rl=1, w_min= -8, w_max=8, w_num=500, dB=True, dB_limit= -100)
#    
#    M = M.real
#    M_orig = M
#    
#    loop = True
#    while(loop):
#        p=input ('Enter the element position p value:  ')
#        q=input('Enter the element position  q value:  ')
#        i=input('Enter the pivot i value:  ')
#        j=input('Enter the pivot j value: ')    
#        
#        tht=annihilation(p,q,i,j,N,M)
#        R=rotation_R(N+2,i,j,tht)
#        
#        T=np.dot(R,M)
#        M=np.dot(T,np.transpose(R))
#        print M        
#        MN2_to_Sparam(M, Rs=1, Rl=1, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)
#        
#        goback = input('If you want to goback to the original M, give 1')
#        if(goback):
#            M = M_orig
#            print M
#        loop = input('If you want to annihilate further, give 1')
#    fbw = input('Specify the value of fractional bandwidth (not pecentage): ')
#    M_final = M*fbw; print 'M_final:','\n', M_final
#    Qi_final = 1/(Rs*fbw); print 'Qi_final:','\n', Qi_final
#    Qo_final = 1/(Rl*fbw); print 'Qo_final:','\n', Qo_final    
    

#==============================================================================
# 3rd order example (zinka)
#==============================================================================

#    N = 4
##    poles = np.array([1.2576, -0.1546 - 0.9218j, -0.1546 + 0.9218j])
#    poles = np.array([-3])
#    eps = 6.0251j; eps_R = 1
##    poles = np.array([])
#    F, P = Chebyshev_gen(N, poles)
##    plot_rational(F, P, x_min= -1, x_max=1, x_num=1000)    
#    F = w_to_s(F, coef_norm=True)
#    P = w_to_s(P, coef_norm=True)      
#    print 'F:', '\n', F; print 'P:', '\n', P
#    [E, roots_E] = poly_E(eps, eps_R, F, P)
#    print 'E:', '\n', E
#    print roots_E # polynomial E has multiple roots ... when N-nfz is odd
#    plot_mag(eps, eps_R, F, P, E, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)
##    plot_delay(roots_E)
#    
#    # From now onwards, unlike the Cameron's example, this filter is doubly terminated
#    M, Rs, Rl = coupling_N(F, P, E, eps, eps_R)
#    print 'M:', '\n', M
#    print 'Rs:', Rs
#    print 'Rl:', Rl
#   
#    MN_to_Sparam(M, Rs, Rl, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)

#==============================================================================
# 3rd order example (P. 336, Sec. 10.3.3, J. S. Hong)
#==============================================================================

#    M = np.array([[-0.27644,1.07534,-0.65769],[1.07534,0.48564,1.07534],[-0.65769,1.07534,-0.27644]])
#    print M
#    MN2_to_Sparam(M, Rs=1/0.69484, Rl=1/0.69484, w_min= -8, w_max=8, w_num=500, dB=True, dB_limit= -50)

#==============================================================================
# 4th order example (IEEE MTT, J.S.Hong, 48, 7, July 2000)
#==============================================================================

#    N = 6
#    poles = np.array([-1.2, 1.2])
#    eps = 2.2j; eps_R = 1
##    poles = np.array([])
#    F, P = Chebyshev_gen(N, poles)
##    plot_rational(F, P, x_min= -1, x_max=1, x_num=1000)
#    F = w_to_s(F, coef_norm=True)
#    P = w_to_s(P, coef_norm=True)    
#    print 'F:', '\n', F; print 'P:', '\n', P
#    print 'E:', '\n', poly_E(eps, eps_R, F, P)[0]
#    E = poly_E(eps, eps_R, F, P)[0]
#    plot_mag(eps, eps_R, F, P, E, w_min= -4, w_max=4)    
#    
##    np.set_printoptions(precision=4, suppress=True)    
##    M = coupling_N2(F, P, E, eps, eps_R)
##    print 'M:', '\n', M.real
##    MN2_to_Sparam(M, Rs=1, Rl=1, w_min= -8, w_max=8, w_num=500, dB=True, dB_limit= -50)
#    
#    # From now onwards, unlike the Cameron's example, this filter is doubly terminated
#    M, Rs, Rl = coupling_N(F, P, E, eps, eps_R)
#    print 'M:', '\n', M
#    print 'Rs:', Rs
#    print 'Rl:', Rl
#   
#    MN_to_Sparam(M, Rs, Rl, w_min= -4, w_max=4, w_num=500, dB=True, dB_limit= -40)
     
#==============================================================================
# Notes to myself
#==============================================================================

# use new "polynomial" class of Numpy in future ...
# always, always, always remember ... eps could be real or imaginary ... depends upon N-nfz

#    N = 6
#    poles = np.array([+1.2])
#    eps = 2.2j; eps_R = 1
##    poles = np.array([])
#    F, P = Chebyshev_gen(N, poles)
##    plot_rational(F, P, x_min= -1, x_max=1, x_num=1000)
#    F = w_to_s(F, coef_norm=True)
#    P = w_to_s(P, coef_norm=True)    
#    print 'F:', '\n', F; print 'P:', '\n', P
#    print 'E:', '\n', poly_E(eps, eps_R, F, P)[0]
#    E = poly_E(eps, eps_R, F, P)[0]
#    plot_mag(eps, eps_R, F, P, E, w_min= -4, w_max=4, dB_limit=-100)    
#    
##    np.set_printoptions(precision=4, suppress=True)    
##    M = coupling_N2(F, P, E, eps, eps_R)
##    print 'M:', '\n', M.real
##    MN2_to_Sparam(M, Rs=1, Rl=1, w_min= -8, w_max=8, w_num=500, dB=True, dB_limit= -50)
#    
#    # From now onwards, unlike the Cameron's example, this filter is doubly terminated
#    M, Rs, Rl = coupling_N(F, P, E, eps, eps_R)
#    print 'M:', '\n', M
#    print 'Rs:', Rs
#    print 'Rl:', Rl
#   
#    MN_to_Sparam(M, Rs, Rl, w_min= -4, w_max=4, w_num=500, dB=True, dB_limit= -100)

#==============================================================================
# SET Conf
#==============================================================================

#    N = 4
##    poles = np.array([1.2576, -0.1546 - 0.9218j, -0.1546 + 0.9218j])
#    poles = np.array([-1.768,1.768])
#    eps = 2.0251j; eps_R = 1
##    poles = np.array([])
#    F, P = Chebyshev_gen(N, poles)
##    plot_rational(F, P, x_min= -1, x_max=1, x_num=1000)    
#    F = w_to_s(F, coef_norm=True)
#    P = w_to_s(P, coef_norm=True)      
#    print 'F:', '\n', F; print 'P:', '\n', P
#    [E, roots_E] = poly_E(eps, eps_R, F, P)
#    print 'E:', '\n', E
#    print roots_E # polynomial E has multiple roots ... when N-nfz is odd
#    plot_mag(eps, eps_R, F, P, E, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)
##    plot_delay(roots_E)
#    
#    # From now onwards, unlike the Cameron's example, this filter is doubly terminated
#    M, Rs, Rl = coupling_N(F, P, E, eps, eps_R)
#    print 'M:', '\n', M.real
#    print 'Rs:', Rs
#    print 'Rl:', Rl
#   
#    MN_to_Sparam(M, Rs, Rl, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)
#    
#
#    M = M.real
#    M_orig = M
#    
#    loop = True
#    while(loop):
#        p=input ('Enter the element position p value:  ')
#        q=input('Enter the element position  q value:  ')
#        i=input('Enter the pivot i value:  ')
#        j=input('Enter the pivot j value: ')    
#        
#        tht=annihilation(p,q,i,j,N,M)
#        R=rotation_R(N,i,j,tht)
#        
#        T=np.dot(R,M)
#        M=np.dot(T,np.transpose(R))
#        print M        
#        MN_to_Sparam(M, Rs, Rl, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)
#        
#        goback = input('If you want to goback to the original M, give 1')
#        if(goback):
#            M = M_orig
#            print M
#        loop = input('If you want to annihilate further, give 1')
#    fbw = input('Specify the value of fractional bandwidth (not pecentage): ')
#    M_final = M*fbw; print 'M_final:','\n', M_final
#    Qi_final = 1/(Rs*fbw); print 'Qi_final:','\n', Qi_final
#    Qo_final = 1/(Rl*fbw); print 'Qo_final:','\n', Qo_final

#    gap, Ke = K(plot=True)
#    gap, Km = K(plot=False)
#    gap, Km2 = K(plot=False)
#    plt.plot(gap, Ke, '-b', label='Electric'); plt.plot(gap, Km,'-r', label='Magnetic')
#    plt.plot(gap, Km2,'-g', label='Magnetic2')
#    plt.axis('tight'); plt.grid(True)
#    plt.xlabel('gap')
#    plt.ylabel('K')
#    plt.legend()
#    plt.show()

#==============================================================================
# Evaluating Coupling Coefficients (3 plots)
#==============================================================================

#    gap, Ke = K(plot=True)
#    gap, Km = K(plot=False)
#    gap, Km2 = K(plot=False)
#    plt.plot(gap, Ke, '-b', label='Electric'); plt.plot(gap, Km,'-r', label='Magnetic')
#    plt.plot(gap, Km2,'-g', label='Magnetic2')
#    plt.axis('tight'); plt.grid(True)
#    plt.xlabel('gap')
#    plt.ylabel('K')
#    plt.legend()
#    plt.show()

#==============================================================================
# Evaluating Coupling Coefficients (1 plot)
#==============================================================================
    
#    gap, K = K(plot=False)
#    plt.plot(gap, K, '-b', label='Magnetic')
#    plt.axis('tight'); plt.grid(True)
#    plt.xlabel('gap')
#    plt.ylabel('K')
#    plt.legend()
#    plt.show()    

#==============================================================================
# Hair-pin filter
#==============================================================================

#     N = 6
# #     poles = np.array([1.2576, -0.1546 - 0.9218j, -0.1546 + 0.9218j])
#     poles = np.array([])
#     eps = 1.6j; eps_R = 1
#      
# #    LAr = 0.5 # pass-band ripple in dB (>0)
# #    LR = -16.426 # return-loss (<0)
# #    print -10*np.log10(1-10**(0.1*LR))
# #    
# #    eps = np.sqrt(10**(LAr*0.1)-1)
# #    print eps
# #    print -10*np.log10(1/(1+eps**2))
# #    print 10*np.log10(1-(1/(1+eps**2)))
#      
# #    poles = np.array([])
#     F, P = Chebyshev_gen(N, poles)   
# #    plot_rational(F, P, x_min= -1, x_max=1, x_num=1000)    
#     F = w_to_s(F, coef_norm=True)
#     P = w_to_s(P, coef_norm=True)      
#     print 'F:', '\n', F; print 'P:', '\n', P
#     [E, roots_E] = poly_E(eps, eps_R, F, P)
#     print 'E:', '\n', E
#     print roots_E # polynomial E has multiple roots ... when N-nfz is odd
#     plot_mag(eps, eps_R, F, P, E, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)
# #    plot_delay(roots_E)
# #    
#     # From now onwards, unlike the Cameron's example, this filter is doubly terminated
#     M, Rs, Rl = coupling_N(F, P, E, eps, eps_R)
#     print 'M:', '\n', M.real
#     print 'Rs:', Rs
#     print 'Rl:', Rl
# #   
# #    MN_to_Sparam(M, Rs, Rl, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)
# #    
# #
#     M = M.real
#     M_orig = M
#      
#     loop = True
#     while(loop):
#         p=input ('Enter the element position p value:  ')
#         q=input('Enter the element position  q value:  ')
#         i=input('Enter the pivot i value:  ')
#         j=input('Enter the pivot j value: ')    
#          
#         tht=annihilation(p,q,i,j,N,M)
#         R=rotation_R(N,i,j,tht)
#          
#         T=np.dot(R,M)
#         M=np.dot(T,np.transpose(R))
#         print M        
#         MN_to_Sparam(M, Rs, Rl, w_min= -10, w_max=10, w_num=500, dB=True, dB_limit= -100)
#          
#         goback = input('If you want to goback to the original M, give 1')
#         if(goback):
#             M = M_orig
#             print M
#         loop = input('If you want to annihilate further, give 1')
#     fbw = input('Specify the value of fractional bandwidth (not pecentage): ')
#     M_final = M*fbw; print 'M_final:','\n', M_final
#     Qi_final = 1/(Rs*fbw); print 'Qi_final:','\n', Qi_final
#     Qo_final = 1/(Rl*fbw); print 'Qo_final:','\n', Qo_final