# -*- coding: utf-8 -*-
"""
Created on Thu Feb 09 23:57:41 2017

@author: Bernt_Lie
"""
import numpy as np
import matplotlib.pyplot as plt
#
NU=np.array([[2130,2300,2357],
    [500,480,600],
    [100,90,150],
    [14.3,11.5,24.3],
    [5.1,3.9,8.],
    [1.9,np.nan,2.6]])
T = np.array([-65,-40,0,100,210,400])
T = (T-32)*5/9;
#
plt.semilogy(T,NU[:,0], 'k-', label="Univis J-43")
plt.semilogy(T,NU[:,1], 'b:', label="Skydrol 500A")
plt.semilogy(T,NU[:,2], 'r--', label="Oronite 8315")
plt.xlabel(r'$T$ [${}^\circ$C]')
plt.ylabel(r'$\nu(T)$ [cSt]')
plt.title(r'Kinematic viscosity $\nu$ as a function of temperature $T$')
plt.grid(True)
plt.legend()