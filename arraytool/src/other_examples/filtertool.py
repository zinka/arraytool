#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

from __future__ import division
import numpy as np
import sympy as sp
import scipy.signal as signal
import matplotlib.pyplot as plt

# adjusting "matplotlib" label fonts ... can't save .svz files using this option
from matplotlib import rc
rc('text', usetex=True)

def cutoff(F, dB_limit= -40):
    r"""
    When magnitude of S11 or S21 is 0, their dB value is '-infinity'. So, this 
    function will be used to cut-off all the value below some 'dB_limit'.
    """
    msk1 = F < dB_limit
    fill = msk1 * dB_limit
    msk2 = F >= dB_limit
    F = F * (msk2) + fill
    return F

def I_to_i(A):
    """Function to convert the symbol 'I' of SymPy to 'i'."""
    B = np.zeros((len(A), 1), dtype=complex)
    for i in range(len(A)):
        B[i] = complex(A[i])
    return B

def s_to_w(poly_ip, coef_norm=False):
    """Function to convert a given polynomial from 's' to 'w' domain (s=jw)."""
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
    """Function to convert a given polynomial from 'w' to 's' domain (w=-js)."""
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
    """Simple function to plot a given rational function (i.e., a ratio of two 
       polynomial functions) of x."""
    x = np.linspace(x_min, x_max, x_num)
    y_Num = np.polyval(Num, x); y_Den = np.polyval(Den, x); y = y_Num / y_Den    
    plt.plot(x, y); plt.axis('tight'); plt.grid(True)
    plt.xlabel('x'); plt.ylabel('y'); plt.show()
    return x, y

def Gramm_Schmidt(A):
    """Function to implement the modified Gram-Schmidt orthonormalization."""
    A = A.T; m, n = A.shape
    q = np.zeros((m, n)); r = np.zeros((n, n))
    for k in range(n):
        r[k, k] = np.linalg.norm(A[:, k])
        q[:, k] = A[:, k] / r[k, k]
        r[k, k + 1:n] = np.dot(q[:, k], A[:, k + 1:n])
        A[:, k + 1:n] = A[:, k + 1:n] - np.outer(q[:, k], r[k, k + 1:n])
    return q.T

def Chebyshev_gen(N, poles):
    """Function to evaluate the numerator and denominator polynomials of the 
       generalized Chebyshev (i.e., with/without poles) filtering function."""
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
    """Function to obtain the polynomial E and its roots in the s-domain."""
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
             dB_limit= -40, show=True):
    """Function to plot magnitudes of S11 and S21 in either linear or dB scale."""
    w = np.linspace(w_min, w_max, w_num)
    F = s_to_w(F); P = s_to_w(P); E = s_to_w(E)
    F_val = np.polyval(F, w); P_val = np.polyval(P, w); E_val = np.polyval(E, w)
    S11 = (1 / eps_R) * (F_val / E_val); S21 = (1 / eps) * (P_val / E_val)
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
    if(show): plt.show()
    return w, S11, S21

def plot_delay(roots_E, w_min= -2, w_max=2, w_num=500, show=True):
    """Function to plot the group delay."""
    w = np.linspace(w_min, w_max, w_num)
    tau = np.zeros_like(w)
    for i in range(len(w)):
        tmp = 0
        for j in range(len(roots_E)):
            tmp += -roots_E[j].real / (roots_E[j].real ** 2 + (w[i] - roots_E[j].imag) ** 2)
        tau[i] = tmp
    plt.plot(w, tau)
    plt.axis('tight'); plt.grid(True)
    plt.xlabel(r'$\Omega\ \mathrm{(rad/s)}$', fontsize=14)
    plt.ylabel(r'$\mathrm{Group}\ \mathrm{delay} \mathrm{\ (s)}$', fontsize=14)    
    if(show): plt.show()        
    return w, tau

def coupling_N(F, P, E, eps, eps_R):
    """Function to evaluate the (N,N) coupling matrix.
       Not done yet... needs to be tested"""
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
    T = np.eye(len(T1k), len(T1k))
    T[0, :] = T1k; T[1, :] = TNk
    np.set_printoptions(precision=6, suppress=True)
    T = Gramm_Schmidt(T) # "normalizing of T1k, TNk" is done in this step
    # swapping the second and last rows after normalization is finished
    temp = np.copy(T[1, :]); T[1, :] = T[-1, :]; T[-1, :] = temp
    Lamb = np.diag(lambdk) # diagonal eigenvalue matrix
    M = np.dot(T, np.dot(Lamb, T.T)) # (N,N) coupling matrix
    return M, RS_L1, RL_LN

def MN_to_Sparam(M, Rs, Rl, w_min= -2, w_max=2, w_num=500, dB=True,
                 dB_limit= -40, plot=True, show=True):
    """Function to plot S parameters from a given (N,N) coupling matrix."""
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
        if(show): plt.show()        
    return w, S11, S21

def coupling_N2(F, P, E, eps, eps_R):
    """Function to evaluate the (N+2,N+2) coupling matrix."""
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
                 dB_limit= -40, plot=True, show=True):
    """Function to plot S parameters from a given (N,N) coupling matrix."""
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
        if(show): plt.show()        
    return w, S11, S21

if __name__ == '__main__':
        
#==============================================================================
# 7th order example (P. 300, Sec. 8.3.1, R. J. Cameron et al.)
#==============================================================================    
    
    N = 7
    poles = np.array([1.2576, -0.1546 - 0.9218j, -0.1546 + 0.9218j])
    eps = 6.0251j; eps_R = 1
#    poles = np.array([])
    F, P = Chebyshev_gen(N, poles)
#    plot_rational(F, P, x_min= -1, x_max=1, x_num=1000)    
    F = w_to_s(F, coef_norm=True)
    P = w_to_s(P, coef_norm=True)      
    print 'F:', '\n', F; print 'P:', '\n', P
    [E, roots_E] = poly_E(eps, eps_R, F, P)
    print 'E:', '\n', E
#    plot_mag(eps, eps_R, F, P, E)
#    plot_delay(roots_E)
    
#    # testing M_to_Sparam function
#    Rs = 1.0442; Rl = 1.0442
#    M = np.array([[0, 0.8577, 0, -0.2174],
#                  [0.8577, 0, 0.7856, 0],
#                  [0, 0.7856, 0, 0.8577],
#                  [-0.2174, 0, 0.8577, 0]]) # a test matrix from a journal paper
#    print 'M:', '\n', M

#    # testing M_to_Sparam function
#    Rs = 1.025; Rl = 1.025
#    M = np.array([[0,-0.863,0,0.02,0],
#                  [-0.863,0,-0.647,0,0],
#                  [0,-0.647,0,-0.632,0],
#                  [0.02,0,-0.632,0,-0.863],
#                  [0,0,0,-0.863,0]]) # another testing coupling matrix
#    print 'M:', '\n', M

    # From now onwards, unlike the Cameron's example, the filter is doubly terminated
    M, Rs, Rl = coupling_N(F, P, E, eps, eps_R)
    print 'M:', '\n', M
    print 'Rs:', Rs
    print 'Rl:', Rl
   
    MN_to_Sparam(M, Rs, Rl, w_min= -2, w_max=2, w_num=500, dB=True, dB_limit= -100)
    
#    # Checking the obtained coupling matrix by plotting the S-parameters
#    M_to_Sparam(M, Rs=1, Rl=1, L1=1, LN=1)

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
  
#==============================================================================
# Notes to myself
#==============================================================================
# use new "polynomial" class of Numpy in future ...