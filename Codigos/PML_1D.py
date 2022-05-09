#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  9 22:01:48 2022

@author: Benjamín Pascual Estrugo, Francisco José Pirt Bastida
"""

import numpy as np
from math import exp
from matplotlib import pyplot as plt
from math import sin, pi
import matplotlib.animation as animation

l = 200

# Parametros del pulso
lc = int(l / 2)
t0 = 40
spread = 12
#Condiciones de contorno
boundary_low = [0, 0]
boundary_high = [0, 0]

nsteps = 305
#Grid
Ex = np.zeros([nsteps,l])
Hy = np.zeros([nsteps,l])
# Loop FDTD
for time_step in range(0, nsteps):
    # Calculamos Ex
    for k in range(1, l):
        Ex[time_step,k] = Ex[time_step-1,k] +  0.5*(Hy[time_step-1,k - 1] - Hy[time_step-1,k])
        
    # Pulso gaussiano en el medio
    pulse = exp(-0.5 * ((t0 - time_step) / spread) ** 2)
    #pulse = 0.5*sin(2 * pi *0.1 * (time_step) / spread)
    Ex[time_step,lc] = pulse
    
    # Condiciones de contorno
    Ex[time_step,0] = boundary_low.pop(0)
    boundary_low.append(Ex[time_step,1])
    
    Ex[time_step,l - 1] = boundary_high.pop(0)
    boundary_high.append(Ex[time_step,l - 2])
    
    # Calculamos Hy
    for k in range(l - 1):
        Hy[time_step,k] = Hy[time_step-1,k] + 0.5 * (Ex[time_step,k] - Ex[time_step,k + 1])
    
           


fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)

def plot_e_field(data, timestep, label):
    ax.plot(data, color='b', linewidth=1)
    ax.set_ylabel('E$_x$', fontsize='14')
    ax.set_xticks(np.arange(0, 199, step=20))
    ax.set_xlim(0, 199)
    ax.set_yticks(np.arange(0, 1.2, step=1))
    ax.set_ylim(-0.2, 1.2)
    ax.set_xlabel('{}'.format(label))
    
    
        
def actualizar(i):
    ax.clear()
    plot_e_field(Ex[2*i],2*i,'celdas')
    

ani = animation.FuncAnimation(fig,actualizar,range(nsteps))
plt.show()



