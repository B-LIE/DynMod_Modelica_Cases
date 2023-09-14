# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 11:19:20 2017

@author: Bernt_Lie
"""
import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sy
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
#s0 = sy.mpmath.taylor(sy.sin,np.pi/2,0)[::-1]
#s1 = sy.mpmath.taylor(sy.sin,np.pi/2,1)[::-1]
#s3 = sy.mpmath.taylor(sy.sin,np.pi/2,3)[::-1]
#s5 = sy.mpmath.taylor(sy.sin,np.pi/2,5)[::-1]
x = sy.Symbol("x")

f = sy.lambdify(x,sy.sin(x),"numpy")
# For n=2, we need a fix for function "lambdify", which 
# doesn't correctly handle constant value functions
e_sin1 = sy.series(sy.sin(x),x0=sy.pi/2,n=2).removeO()
f_sin1px = sy.lambdify(x,e_sin1+x,"numpy")
f1 = lambda x: f_sin1px(x)-x
# Remaining cases work with "lambdify" function                      
f2 = sy.lambdify(x,sy.series(sy.sin(x),x0=sy.pi/2,n=3).removeO(),"numpy")
f3 = sy.lambdify(x,sy.series(sy.sin(x),x0=sy.pi/2,n=4).removeO(),"numpy")
f4 = sy.lambdify(x,sy.series(sy.sin(x),x0=sy.pi/2,n=5).removeO(),"numpy")
f5 = sy.lambdify(x,sy.series(sy.sin(x),x0=sy.pi/2,n=6).removeO(),"numpy")
#
x=np.linspace(0,2*np.pi)
#
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,f1(x), color=Cb1,ls=LS3, linewidth=LW1, label="0,1 order")
#plt.plot(x,f2(x), 'g--D', linewidth=2, label="2 order")
ax.plot(x,f3(x), color=Cg1,ls=LS3, linewidth=LW1, label="2,3 order")
#plt.plot(x,f4(x), 'k-^', linewidth=2, label="4 order")
ax.plot(x,f5(x), color=Cr1,ls=LS3, linewidth=LW1, label="4,5 order")
ax.plot(x,np.sin(x),'k-', linewidth=LW1, label=r"$\sin(x)$")
ax.set_xlabel(r"$x$")
ax.set_ylim(-1.5,1.5)
ax.set_xlim(0,2*np.pi)
ax.grid()
ax.legend()
dx = np.pi/4
xticks = [i*dx for i in range(9)]
xtick_labels = [r"$0$", r"$\pi/4$"] + [r"${}\pi/4$".format(k) for k in range(2,9)]
ax.set_xticks(xticks)
ax.set_xticklabels(xtick_labels)
#
fig_name = "sin-approx.pdf"
plt.savefig(fig_path+fig_name)