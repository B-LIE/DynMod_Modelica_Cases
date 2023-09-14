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
x = lambda t: 1.-t-100.*np.sin(t/10)+10.*t*np.cos(t/10.)
v = lambda t: 9.-10.*np.cos(t/10.)
#
t=np.linspace(0,5.)
#
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t,x(t), color=Cb1,ls=LS1, linewidth=LW1, label=r"$x$")
ax.plot(t,v(t), 'r--', linewidth=LW1, label=r"$v$")
#
ax.set_xlabel(r"$t$")
ax.set_xlim(0,5.)
ax.set_title(r"Solutions, non-diagonalizable system with input $u(t) = \sin(\frac{t}{10})$")
ax.grid()
ax.legend()
#
fig_name = "sin-input-non-diagonalizable.pdf"
plt.savefig(fig_path+fig_name)