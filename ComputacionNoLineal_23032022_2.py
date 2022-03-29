#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 13:50:35 2022

@author: juanjosefernandezmorales
"""

import numpy as np
import matplotlib.pyplot as plt

# -------------------------- Datos

# Entradas Posiciones y tiempos
xmax = 10
xmin = 0
N    = 101

# Entradas del campo 
mu = 1
ep = 1
speedoflight = 1/np.sqrt(mu*ep)

t = 0
tend = 2.5

# Entradas gaussiana
mean   = (xmax + xmin)/2
spread = (xmax - xmin)/10

# grid y dualgrid
grid     = np.linspace(xmin, xmax, N)                                          
dualgrid = (grid[:-1] + grid[1:])/2                                            

dx = grid[1] - grid[0]
dt = 0.8*dx/speedoflight


# ---------------------------- Calculo de los campos

E = np.zeros((round(tend/dt), N))
H = np.zeros((round(tend/dt), N-1))

E[0] = np.exp(- np.power(grid - mean, 2) / 2/ np.sqrt(spread))     
H[0] = dualgrid * 0

for i in range(1, (round(tend/dt))):
    E[i][1:-1] = - (dt/dx) * (H[i-1][1:] - H[i-1][:-1] + E[i-1][1:-1])
    H[i] = - (dt/dx) * (E[i-1][1:] - E[i-1][:-1] + H[i-1])


# ---------------------------- Animacion

import matplotlib.animation as animation

fig, ax = plt.subplots()

xEdata, yEdata = [], []
lnE, = plt.plot([], [], '-')

xHdata, yHdata = [], []
lnH, = plt.plot([], [], '-')

def init():
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-1, 1)
    return lnE, lnH

def update(frame):
    xEdata = grid
    yEdata = E[frame]
    lnE.set_data(xEdata, yEdata)
    
    xHdata = grid[1:]
    yHdata = H[frame]
    lnH.set_data(xHdata, yHdata)
    
    return lnE, lnH

ani = animation.FuncAnimation(fig, update, frames=range((round(tend/dt))),
                    init_func=init, blit=True)


ani.save(r"prueba.gif")