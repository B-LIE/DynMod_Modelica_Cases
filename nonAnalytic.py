# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 10:22:44 2017

@author: Bernt_Lie
"""
import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt
import pandas as pd
LW = 2
#
def non_analytic(x):
    val = np.zeros(x.shape)
    ineg = x<0
    val[np.nonzero(ineg)] = 0
    val[np.nonzero(~ineg)] = x[np.nonzero(~ineg)]**2
    return val

x = np.linspace(-5,5)
plt.plot(x,non_analytic(x), linewidth=LW)
plt.xlabel(r"$x$")
plt.ylabel(r"$f(x)$")
plt.title(r"Example of non-analytic function")