#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 19:34:53 2022

@author: benjapases (Benjamín Pascual Estrugo)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#----------------
xMax = 10.0
xMin=0.0
N = 101

mean = (xMax + xMin)/2#media de la gaussiana
spread = ((xMax - xMin)/10)

tEnd = 2.5



speedOfLight = 1.0

#----------------

grid = np.linspace(xMin,xMax,N)#grid primario
dualGrid = (grid[:-1]+grid[1:])/2#grid dual
dx = grid[1]-grid[0]#cualquiera dos puntos consecutivos serviría
dt = 0.8*dx/speedOfLight
dimt = round(tEnd/dt)#numero de filas (tiempos)

Eold = np.exp( -np.power(grid - mean,2)/2.0/np.sqrt((spread)))#condición inicial
Hold = dualGrid*0.0#condición lineal para H

Enew = np.zeros([dimt,101])
Hnew = np.zeros([dimt,100])

Enew[0][:] = Eold



t = 0.0
i = 1
while t<tEnd:

 Enew[i][1:-1] = -(dt/dx)*(Hnew[i-1][1:] - Hnew[i-1][:-1]) + Enew[i-1][1:-1]
 Hnew[i][:] = - (dt/dx)*(Enew[i-1][1:] - Enew[i-1][:-1]) + Hnew[i-1][:]

 t += dt
 i += 1
 if i == 31:
     break


#Create figure
fig = plt.figure(figsize = (15,15))

ax = fig.add_subplot(111)

#funcion animacion
def actualizar(i):
    ax.clear()
    ax.plot(grid,Enew[i,:],color = "darkblue")
    ax.plot(dualGrid,Hnew[i,:], color = "red")
    ax.set_ylim(-1,1)#esto hace que el eje y no se mueva en la animación




ani = animation.FuncAnimation(fig,actualizar,range(dimt))




















