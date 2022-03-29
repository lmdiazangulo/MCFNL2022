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

tEnd = 10.0



speedOfLight = 1.0

#----------------

grid = np.linspace(xMin,xMax,N)#grid primario
dualGrid = (grid[:-1]+grid[1:])/2#grid dual
dx = grid[1]-grid[0]#cualquiera dos puntos consecutivos serviría
dt = 0.8*dx/speedOfLight
dimt = round(tEnd/dt)+2#numero de filas (tiempos)

Eold = np.exp( -np.power(grid - mean,2)/2.0/np.sqrt((spread)))#condición inicial
Hold = dualGrid*0.0#condición lineal para H

E = np.zeros([dimt,101])
H = np.zeros([dimt,100])

E[0][:] = Eold



t = 0.0
i = 1
while t<tEnd:

    E[i][1:-1] = -(dt/dx)*(H[i-1][1:] - H[i-1][:-1]) + E[i-1][1:-1]
    E[i][-1] = -(dt/dx)*(- 2.0 * H[i-1][-1])+  E[i-1][-1]
    H[i][:] = - (dt/dx)*(E[i][1:] - E[i][:-1]) + H[i-1][:]

    t += dt
    i += 1


#Create figure
fig = plt.figure(figsize = (15,15))

ax = fig.add_subplot(111)

#funcion animacion
def actualizar(i):
    ax.clear()
    ax.plot(grid,E[i,:],color = "darkblue")
    ax.plot(dualGrid,H[i,:], color = "red")
    ax.set_ylim(-1,1)#esto hace que el eje y no se mueva en la animación




ani = animation.FuncAnimation(fig,actualizar,range(dimt))



plt.show()
print("END")
















