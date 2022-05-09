#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot    as plt
import matplotlib.animation as animation
from fdtd   import FDTD, Source
from panel  import RTfunction

#----------------
# Grid Properties
N     = 200
dx    = 0.002
xMax  = N*dx
tEnd  = 5e-8

# Panel Properties
sigma        = 0.04
permittivity = 4
panel_min    = 110
panel_max    = 140
thickness    = (panel_min - panel_max)*dx

# Pulse Properties


#----------------

grid  = np.linspace(0, dx*N, N)

fdtd1 = FDTD(grid_prop  = (dx, N, tEnd), 
             panel_prop = (panel_min, panel_max, sigma, permittivity),
             pulse  = Source('gauss', 40, 12, 20))

fdtd2 = FDTD(grid_prop  = (dx, N, tEnd), 
             panel_prop = (panel_min, panel_max, 0, 1),
             pulse  = Source('gauss', 40, 12, 20))

probedE1, probedH1 = fdtd1.simulation()
probedE2, probedH2 = fdtd2.simulation() 

R, T, w = RTfunction(probedE1, probedE2, panel_min, panel_max, tEnd)     

#Create figure
fig, ax = plt.subplots(figsize=(10, 5))
ax.set(xlim=(0, N), ylim=(-1.2, 1.2))

#funcion escenario
def init():
    pass

def panel():
    plt.hlines(0, 0, panel_min, 'k')
    plt.hlines(0, panel_max, N, 'k')
    plt.vlines(panel_min, 0, 1, 'k')
    plt.vlines(panel_max, 0, 1, 'k')

#funcion animacion
def actualizar(i):
    ax.clear()
    ax.plot(np.arange(0, N), probedE1[i,:], ".", markersize=5)
    ax.plot(np.arange(0, N), probedH1[i,:], ".", markersize=5)
    panel()
    ax.grid()
    ax.set_ylim(-1,1)

# Panel
ani = animation.FuncAnimation(fig, actualizar,
                              frames = range(probedE1.shape[0]), interval= 1)
plt.draw()
plt.show()

plt.cla()

plt.figure()
plt.xlim((1,7e10))
plt.ylim((0,1))

plt.plot(w, R, label='FDTD R')
plt.grid()
plt.legend(loc='best', ncol=4)
plt.show()
plt.show()

print("END")















