# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 09:45:30 2017

@author: Bernt_Lie
"""
import numpy as np
import matplotlib.pyplot as plt
#
X = np.linspace(1e-3,1e-2)
Y1 = np.zeros(np.shape(X))
idx = range(0,len(Y1))
for i in idx:
    if X[i] < 0.005:
        Y1[i] = 0.7*0.024*X[i]
    else:
        Y1[i] = 0.7*0.024*(0.005-(X[i]-0.005))
Y2 = 0.6*0.024*X/np.sqrt(1-(0.6*X/0.01)**2)
#
plt.plot(X,Y1,label="Kurode et al.")
plt.plot(X,Y2,label="Merritt")

plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$C_\mathrm{d}A$')
plt.title(r'Valve coefficient model')
plt.grid(True)
plt.legend()