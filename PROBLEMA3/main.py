#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot    as plt

from solver.fdtd   import FDTD, Source
from solver.panel  import RTfunction, Panel


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
thickness    = (panel_max - panel_min)*dx

###############################################################################
# FDTD
###############################################################################

grid  = np.linspace(0, dx*N, N)

# Simulación principal
fdtd1 = FDTD(grid_prop  = (dx, N, tEnd), 
             panel_prop = (panel_min, panel_max, sigma, permittivity),
             pulse  = Source('gauss', 40, 12, 20))

# Simulación vacío
fdtd2 = FDTD(grid_prop  = (dx, N, tEnd), 
             panel_prop = (panel_min, panel_max, 0, 1),
             pulse  = Source('gauss', 40, 12, 20))

probedE1, probedH1 = fdtd1.simulation()
probedE2, probedH2 = fdtd2.simulation() 

R1, T1, W1 = RTfunction(probedE1, probedE2, panel_min, panel_max, tEnd)

###############################################################################
# Panel
###############################################################################    

w = np.linspace(1, 3e12, 100001) * 2 *np.pi
r = np.abs(Panel(thickness,  permittivity,   sigma).R(w))
t = np.abs(Panel(thickness,  permittivity,   sigma).T(w))


###############################################################################
# Visualization Principal
###############################################################################

fig, axs = plt.subplots(1, 1, constrained_layout=True, figsize=(10,5))

axs.set_title('Reflexion and Transmisition coefficients')
axs.set_ylabel(r"$R$ $T$", fontsize = 16)
axs.set_xlabel(r" $\omega$ $[Hz]$", fontsize = 16)


axs.set_xlim(1, 7e10)
axs.set_ylim(0, 1.2)

axs.plot(W1, R1, label='Numerical R', color = "blue")
axs.plot(w, r, label='Theoretical R', color = "orange")

axs.plot(W1, T1, 'b--', label='Numerical T')
axs.plot(w, t, '--', color = "orange", label='Theoretical T')

axs.grid()
axs.legend(loc='best', ncol=4, fontsize=14)

plt.show()

print("END")
