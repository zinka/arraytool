#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2011 Srinivasa Rao Zinka
# License: New BSD License.

"""
In the previous section, plotting in :math:`(u,v)` is described. However, some times
it may be necessary to plot in :math:`(\\theta,\phi)` domain too. In such cases,
one can use the function :func:`pattern_tp`. A simple example to plot
radiation pattern of a planar uniform array in :math:`(\\theta,\phi)` domain is 
shown below::

    import arraytool.planar as planar
    import numpy as np
    
    a = 0.6 # separation between the elements along x-axis (normalized WRS wavelength)
    b = 0.5 # separation between the elements along y-axis (normalized WRS wavelength)
    M = 10 # no. of elements along x-axis
    N = 11 # no. of elements along y-axis
    
    A = np.ones((N, M)) # Uniform planar excitation
    
    # Converting the 'excitation & position' information into Arraytool input format
    array_ip = planar.ip_format(a, b, A)
    
    # Calling the 'pattern_tp' function to evaluate and plot 3D AF/GF/NF    
    [tht, phi, AF] = planar.pattern_tp(array_ip, tht_scan=(0)*np.pi, phi_scan=(0)*np.pi, tht_min= 0,
                      tht_max=0.5*np.pi, tht_num=200, phi_min= 0*np.pi, phi_max=2*np.pi,
                      phi_num=200, scale="dB", dB_limit= -40, factor="GF", plot_type="polar")
                      
.. image:: _static/pattern_tp_1.png
   :align: center                      
                      
If a rectangular surf plot is needed, then the following command can be used::

    # Calling the 'pattern_tp' function to evaluate and plot 3D AF/GF/NF    
    [tht, phi, AF] = planar.pattern_tp(array_ip, tht_scan=(0)*np.pi, phi_scan=(0)*np.pi, tht_min= 0,
                      tht_max=0.5*np.pi, tht_num=200, phi_min= 0*np.pi, phi_max=2*np.pi,
                      phi_num=200, scale="dB", dB_limit= -40, factor="GF", plot_type="rect")

.. image:: _static/pattern_tp_2.png
   :align: center   
               
If a contour plot is needed, then the following command can be used::

    # Calling the 'pattern_tp' function to evaluate and plot 3D AF/GF/NF    
    [tht, phi, AF] = planar.pattern_tp(array_ip, tht_scan=(0)*np.pi, phi_scan=(0)*np.pi, tht_min= 0,
                      tht_max=0.5*np.pi, tht_num=200, phi_min= 0*np.pi, phi_max=2*np.pi,
                      phi_num=200, scale="dB", dB_limit= -40, factor="GF", plot_type="contour")
                      
.. image:: _static/pattern_tp_3.png
   :align: center   
"""

#import arraytool.planar as planar
#import numpy as np
#
#a = 0.6 # separation between the elements along x-axis (normalized WRS wavelength)
#b = 0.5 # separation between the elements along y-axis (normalized WRS wavelength)
#M = 10 # no. of elements along x-axis
#N = 11 # no. of elements along y-axis
#
#A = np.ones((N, M)) # Uniform planar excitation
#
## Converting the 'excitation & position' information into Arraytool input format
#array_ip = planar.ip_format(a, b, A)
#
## Calling the 'pattern_tp' function to evaluate and plot 3D AF/GF/NF    
#[tht, phi, AF] = planar.pattern_tp(array_ip, tht_scan=(0)*np.pi, phi_scan=(0)*np.pi, tht_min= 0,
#                  tht_max=0.5*np.pi, tht_num=200, phi_min= 0*np.pi, phi_max=2*np.pi,
#                  phi_num=200, scale="dB", dB_limit= -40, factor="GF", plot_type="polar")
#
### Calling the 'pattern_tp' function to evaluate and plot 3D AF/GF/NF    
##[tht, phi, AF] = planar.pattern_tp(array_ip, tht_scan=(0)*np.pi, phi_scan=(0)*np.pi, tht_min= 0,
##                  tht_max=0.5*np.pi, tht_num=200, phi_min= 0*np.pi, phi_max=2*np.pi,
##                  phi_num=200, scale="dB", dB_limit= -40, factor="GF", plot_type="rect")
#
### Calling the 'pattern_tp' function to evaluate and plot 3D AF/GF/NF    
##[tht, phi, AF] = planar.pattern_tp(array_ip, tht_scan=(0)*np.pi, phi_scan=(0)*np.pi, tht_min= 0,
##                  tht_max=0.5*np.pi, tht_num=200, phi_min= 0*np.pi, phi_max=2*np.pi,
##                  phi_num=200, scale="dB", dB_limit= -40, factor="GF", plot_type="contour")