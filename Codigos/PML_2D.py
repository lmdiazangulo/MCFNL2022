#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  9 22:13:08 2022

@author: Benjamín Pascual Estrugo, Francisco José Pirt Bastida
"""

import numpy as np
from math import sin, pi
from matplotlib import pyplot as plt
from matplotlib import cm
import mpl_toolkits.mplot3d.axes3d
import matplotlib.animation as animation

#Definimos nuestra malla
ie = 60
je = 60
ic = int(ie / 2 - 5)
jc = int(je / 2 - 5)

ez = np.zeros((ie, je))
dz = np.zeros((ie, je))
hx = np.zeros((ie, je))
hy = np.zeros((ie, je))
ihx = np.zeros((ie, je))
ihy = np.zeros((ie, je))

ddx = 0.01 # Tamaño de la celda
dt = ddx / 6e8 # Time step size

# Create Dielectric Profile
epsz = 8.854e-12
# Parametros del pulso
t0 = 40
spread = 12

gaz = np.ones((ie, je))

# Calculamos los parametros del PML
gi2 = np.ones(ie)
gi3 = np.ones(ie)
fi1 = np.zeros(ie)
fi2 = np.ones(ie)
fi3 = np.ones(ie)

gj2 = np.ones(ie)
gj3 = np.ones(ie)
fj1 = np.zeros(ie)
fj2 = np.ones(ie)
fj3 = np.ones(ie)

# PML
npml = 6
for n in range(npml):
    xnum = npml - n
    xd = npml
    xxn = xnum / xd
    xn = 0.33 * xxn ** 3
    
    gi2[n] = 1 / (1 + xn)
    gi2[ie - 1 - n] = 1 / (1 + xn)
    gi3[n] = (1 - xn) / (1 + xn)
    gi3[ie - 1 - n] = (1 - xn) / (1 + xn)
    
    gj2[n] = 1 / (1 + xn)
    gj2[je - 1 - n] = 1 / (1 + xn)
    gj3[n] = (1 - xn) / (1 + xn)
    gj3[je - 1 - n] = (1 - xn) / (1 + xn)
    
    xxn = (xnum - 0.5) / xd
    xn = 0.33 * xxn ** 3
    
    fi1[n] = xn
    fi1[ie - 2 - n] = xn
    fi2[n] = 1 / (1 + xn)
    fi2[ie - 2 - n] = 1 / (1 + xn)
    fi3[n] = (1 - xn) / (1 + xn)
    fi3[ie - 2 - n] = (1 - xn) / (1 + xn)
    
    fj1[n] = xn
    fj1[je - 2 - n] = xn
    fj2[n] = 1 / (1 + xn)
    fj2[je - 2 - n] = 1 / (1 + xn)
    fj3[n] = (1 - xn) / (1 + xn)
    fj3[je - 2 - n] = (1 - xn) / (1 + xn)
    
nsteps = 200


ez = np.zeros([nsteps,ie,je])

#FDTD Loop
for time_step in range(0, nsteps):
    # Calculamos Dz
    for j in range(1, je):
        for i in range(1, ie):
            dz[i, j] = gi3[i] * gj3[j] * dz[i, j] + \
                       gi2[i] * gj2[j] * 0.5 * \
                           (hy[i, j] - hy[i - 1, j] -
                            hx[i, j] + hx[i, j - 1])
                           
                           
    # Ponemos el pulso justo en medio
    pulse = sin(2 * pi * 1500 * 1e6 * dt * time_step)
    #pulse = np.exp(-0.5*((t0-time_step)/spread)**2)
    dz[ic, jc] = pulse
    
    ez[time_step] = gaz * dz # Calculamos Ez
    
    # Calculamos Hx
    for j in range(je - 1):
        for i in range(ie - 1):
            curl_e = ez[time_step,i, j] - ez[time_step,i, j + 1]
            ihx[i, j] = ihx[i, j] + curl_e
            hx[i, j] = fj3[j] * hx[i, j] + fj2[j] * \
                        (0.5 * curl_e + fi1[i] * ihx[i, j])      
                        
                        
                        
                        
    # Calculamos Hy
    for j in range(0, je - 1):
        for i in range(0, ie - 1):
            curl_e = ez[time_step,i, j] - ez[time_step,i + 1, j]
            ihy[i, j] = ihy[i, j] + curl_e
            hy[i, j] = fi3[i] * hy[i, j] - fi2[i] * \
                (0.5 * curl_e + fj1[j] * ihy[i, j])
                
                
   
plt.rcParams['font.size'] = 12
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.linestyle'] = 'dotted'
fig = plt.figure(figsize=(15, 15))
ax=fig.add_subplot(111,projection="3d")

X, Y = np.meshgrid(range(je), range(ie))


def plot_e_field(ax, data, timestep):
    ax.set_zlim(-0.5, 0.5)
    ax.view_init(elev=30., azim=-135)
    surf = ax.plot_surface(X, Y, data, rstride=1, cstride=1,
                    cmap=cm.coolwarm,
                    edgecolor='black', linewidth=.25)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(r' $E_{Z}$', rotation=90, labelpad=10,
                  fontsize=14)
    ax.set_xlabel('celdas')
    ax.set_ylabel('celdas')
    ax.set_xticks(np.arange(0, 60, step=20))
    ax.set_yticks(np.arange(0, 60, step=20))
    ax.set_zticks([-0.5, 0, 0.5])
    ax.text2D(0.6, 0.7, "T = {}".format(timestep),
              transform=ax.transAxes)
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = \
        ax.zaxis.pane.fill = False
    plt.gca().patch.set_facecolor('white')
    ax.dist = 11
    
   
def actualizar(i):
    ax.clear()
    plot_e_field(ax,ez[i],i)
    
    
ani = animation.FuncAnimation(fig,actualizar,range(200))
  
    
