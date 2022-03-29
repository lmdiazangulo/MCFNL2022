#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 13:50:35 2022

@author: juanjosefernandezmorales
"""

import numpy as np
import matplotlib.pyplot as plt

# --------------------------

# Entradas Posiciones y tiempos
xmax = 10
xmin = 0
N    = 101

# Entradas del campo 
mu = 1
ep = 1
speedoflight = 1/np.sqrt(mu*ep)

t = 0
tend = 0.2

# Entradas gaussiana
mean = (xmax + xmin)/2
spread = (xmax - xmin)/10

# --------------------------

grid     = np.linspace(xmin, xmax, N)                                          # Valores de las posiciones
dualgrid = (grid[:-1] + grid[1:])/2                                            # Valores intermedios posiciones

dx = grid[1] - grid[0]
dt = 0.1*dx/speedoflight


# ----------------------------

Eold = np.exp(- np.power(grid - mean, 2) / 2/ np.sqrt(spread))                 # Gaussiana inicial
Hold = dualgrid * 0

"""
# Visualización de la Gaussiana
plt.plot(grid, Eold)
plt.show()
"""

Enew = Eold * 0
Hnew = Hold * 0

while (t < tend):
    # Existe condiciones de frontera E[0] = E[-1] = 0
    Enew[1:-1] = - (dt/dx) * (Hold[1:] - Hold[:-1] + Eold[1:-1])
    Hnew = - (dt/dx) * (Eold[1:] - Eold[:-1] + Hold)
    
    Eold = Enew
    Eold[0]  = 0
    Eold[-1] = 0
    Hold = Hnew
    t += dt
    
plt.plot(grid, Enew)
plt.plot(grid[1:], Hnew)
plt.show()
plt.cla()

import matplotlib.animation as animation

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'ro')

def init():
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(-1, 1)
    return ln,

def update(frame):
    xdata = [frame]
    ydata = [np.sin(frame)]
    ln.set_data(xdata, ydata)
    return ln,

ani = animation.FuncAnimation(fig, update, frames=np.linspace(0, 2*np.pi, 128),
                    init_func=init, blit=True)


ani.save(r"prueba.gif" )