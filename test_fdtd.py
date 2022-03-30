#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 19:34:53 2022

@author: benjapases (Benjamín Pascual Estrugo)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from solver.fdtd import *

#----------------
xMax = 10.0
xMin=0.0
N = 101

mean = (xMax + xMin)*1/3
spread = ((xMax - xMin)/10)

tEnd = 10.0

#----------------

grid = np.linspace(xMin,xMax,N)
fdtd = FDTD(grid)

fdtd.E = np.exp( -np.power(fdtd.grid - mean,2)/2.0/np.sqrt((spread)))

t = 0.0
probedE = fdtd.E
probedH = fdtd.H
while t<tEnd:
    t = fdtd.step(t)
    probedE = np.vstack((probedE, fdtd.E))
    probedH = np.vstack((probedH, fdtd.H))
    


#Create figure
fig = plt.figure(figsize = (15,15))

ax = fig.add_subplot(111)

#funcion animacion
def actualizar(i):
    ax.clear()
    ax.plot(fdtd.grid,     probedE[i,:], ".")
    ax.plot(fdtd.dualGrid, probedH[i,:], ".")
    ax.grid()
    ax.set_ylim(-1,1)

ani = animation.FuncAnimation(fig,actualizar,range(probedE.shape[0]))

plt.show()
print("END")
















