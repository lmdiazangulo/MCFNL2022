#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from solver.fdtd_nonuniform import *


#----------------
#define the grid
xMax = 10.0
xMin=0.0

#create a non-uniform grid
x = xMin
v = [x]
step = 0.05

while x < xMax:
    r = round(random.uniform(0.0, 1.0), 2)
    x = x+step
    v.append(x)
    step = step/0.95
    
#a = np.linspace(xMin, xMax/3, 10, endpoint=False)
#b = np.linspace(xMax/3, 2*xMax/3, 30, endpoint=False)
#c = np.linspace(2*xMax/3, xMax, 10)
#v = np.concatenate((a,b,c))  

grid = np.array(v)
# use geomspace or logspace
fdtd = FDTD(grid)
mean = (xMax + xMin)*1/3
spread = ((xMax - xMin)/10)
fdtd.E = np.exp( -np.power(fdtd.grid 
   - mean,2)/2.0/np.sqrt((spread)))

t = 0.0
tEnd = 20.0

probedE = fdtd.E
probedH = fdtd.H
while t<tEnd:
    t = fdtd.step(t)
    probedE = np.vstack((probedE, fdtd.E))
    probedH = np.vstack((probedH, fdtd.H))
    
    
#create figure
fig = plt.figure(figsize = (7,5))
ax = fig.add_subplot(111)

#update
def actualizar(i):
    ax.clear()
    ax.plot(fdtd.grid,     probedE[i,:], ".", color='teal')
    ax.plot(fdtd.dualGrid, probedH[i,:], ".", color='violet')
    ax.grid()
    ax.set_ylim(-1,1.2)
    
ani = animation.FuncAnimation(fig,actualizar,range(probedE.shape[0]))
plt.show()
print("END")    
