# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 10:18:45 2017

@author: Bernt_Lie
"""
import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt
import pandas as pd
plt.rc('text', usetex=True)
LW1 = 2.5
LW2 = LW1/2
Cb1 = (0.3,0.3,1)
Cb2 = (0.7,0.7,1)
Cg1 = (0,0.6,0)
Cg2 = (0.5,0.8,0.5)
Cr1 = "Red"
Cr2 = (1,0.5,0.5)
LS1 = "solid"
LS2 = "dotted"
LS3 = "-."
fig_path = "C:/Users/Bernt_Lie/OneDrive/Documents/booksBLSOL/LyX-test/figs/"
#
x1 = lambda t: np.exp(-10*t) + np.exp(10*t)*(1-30*np.exp(-11*t)+29*np.exp(-10*t))/40
x2 = lambda t: -np.exp(-10*t) + np.exp(10*t)*(-1+10*np.exp(-11*t)-11*np.exp(-10*t))/40
#
t=np.linspace(0,0.6)
#
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t,x1(t), color=Cb1,ls=LS1, linewidth=LW1, label=r"$x_1$")
ax.plot(t,x2(t), 'r--', linewidth=LW1, label=r"$x_2$")
#
ax.set_xlabel(r"$t$")
ax.set_xlim(0,0.6)
ax.set_title(r"Solutions with constant input $u(t) = 1$")
ax.grid()
ax.legend()
#
fig_name = "Constant-input.pdf"
plt.savefig(fig_path+fig_name)