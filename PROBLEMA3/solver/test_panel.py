#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from panel import Panel

#Data input
N    = 101
dx   = 0.01
xMax = N*dx
tEnd = 2e-8

sigma        = 0.04
permittivity = 4
panel_min    = 80
panel_max    = 85
thickness    = (panel_max - panel_min)*dx

#Analytical Result
omega = np.linspace(1, 3e12, 100001) * 2 *np.pi
R = np.abs(Panel(thickness,  permittivity,   sigma).R(omega))
T = np.abs(Panel(thickness,  permittivity,   sigma).T(omega))
    
#Plot results
plt.figure()
plt.xlim((1,7e10))
plt.ylim((0,1))

plt.plot(omega, R, label='Exact R')
plt.plot(omega, T, label='Exact T')

plt.grid()
plt.legend(loc='best', ncol=4)
plt.show()
