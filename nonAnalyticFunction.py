# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 10:22:44 2017

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
fig_path = "C:/Users/Bernt_Lie/OneDrive/Documents/booksBLSOL/LyX-test/figs/"
#
def non_analytic(x):
    val = np.zeros(x.shape)
    ineg = x<0
    val[np.nonzero(ineg)] = 0
    val[np.nonzero(~ineg)] = x[np.nonzero(~ineg)]**2
    return val

x = np.linspace(-5,5)
plt.plot(x,non_analytic(x), color=Cb1,linewidth=LW1)
plt.xlabel(r"$x$")
plt.ylabel(r"$f(x)$")
plt.title(r"Non-analytic function")
plt.grid()
fig_name = "nonAnalyticFunction.pdf"
plt.savefig(fig_path+fig_name)