# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 16:02:20 2017

@author: Bernt Lie
"""
import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt
import pandas as pd
#
rv = 1.5
Rv = 2
rY = 1.6
RY = 3
#
Y0 = np.sqrt(RY**2-rY**2)
theta0 = np.arccos(rY/RY)
d0 = Rv-rv
#
L = np.sqrt(d0**2/2)*1.8
alpha10 = np.arccos(d0/(2*L))
#
Y = np.linspace(2.25,3.22)
#
theta = np.arccos((rY**2 + RY**2 - Y**2)/2/rY/RY)
dtheta = theta-theta0
d = np.sqrt(rv**2 + Rv**2-2*rv*Rv*np.cos(dtheta))
psi = np.arccos((d**2 + Rv**2 - rv**2)/(2*d*Rv))*np.sign(dtheta)
phi = np.arccos(d/2/L)
#
alpha1 = phi - psi
#
plt.plot(Y,alpha1*180/np.pi,linewidth=2.5)
plt.xlabel(r'$Y$ [m]')
plt.ylabel(r'$\alpha_1$ [${}^\circ$]')
plt.title(r'Guide Vane angle $\alpha_1$ as a function of actuator position $Y$')
plt.grid(True)
plt.hold(True)
plt.plot(Y0,alpha10*180/np.pi,'ko')
plt.hold(False)