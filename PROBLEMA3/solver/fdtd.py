#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import copy
import numpy as np
import scipy.constants as sp


class FDTD:
    def __init__(self, grid_prop, panel_prop, pulse):
        self.dx, self.ncells, self.time = grid_prop
        self.dt = self.dx / (2*sp.c)
        self.nsteps= int(self.time  / self.dt)
        
        self.pulse = pulse
        
        self.E=np.zeros(self.ncells)
        self.H=np.zeros(self.ncells)

        self.grid_E = np.ones(self.ncells)
        self.grid_H = np.ones(self.ncells) / 2.0
        
        # Panel
        self.panel_min, self.panel_max, self.sigma, self.epsilon_r = panel_prop
        eaf = self.dt * self.sigma / (2 * sp.epsilon_0 * self.epsilon_r)
        
        self.grid_E[self.panel_min: self.panel_max] = (1 - eaf ) / (1 + eaf )
        self.grid_H[self.panel_min: self.panel_max] = 0.5 / (self.epsilon_r * (1 + eaf ))
        
        
        self.ca = self.grid_E[1:-1]
        self.cb = self.grid_H[1:-1]
        
        # Simulation
    def simulation(self):
        probedE = np.zeros_like(self.E)
        probedH = np.zeros_like(self.H)
        
        for time_step in range(1, self.nsteps + 1):
            probedE = np.vstack((probedE, self.E))
            probedH = np.vstack((probedH, self.H))
            self.step(time_step)
            
        probedE = np.vstack((probedE, self.E))
        probedH = np.vstack((probedH, self.H))
    
        
        return probedE, probedH
    
    def step(self, time_step):
        E  = self.E
        H  = self.H
        ca = self.ca
        cb = self.cb
      
        E_old=copy.deepcopy(E)
        
        _ = ca * E[1:-1] + cb * (H[:-2] - H[1:-1])
        E[1:-1] = ca * E[1:-1] + cb * (H[:-2] - H[1:-1])
       
        E[self.pulse.k_ini] +=  0.5*self.pulse.pulse(time_step) 
        
        #boundarymur
        #######################################################################
        ncells, dx, dt= self.ncells, self.dx, self.dt

        c_bound=(sp.c * dt - dx)/(sp.c * dt + dx)

        E[0]=E_old[1] + c_bound * (E[1]-E_old[0])
        E[ncells-1]=E_old[ncells-2] + c_bound * (E[ncells-2]-E_old[ncells-1])
        #######################################################################
        
        H[:-1] = H[:-1] + 0.5 * (E[:-1] - E[1:])   

        t = time_step+1/2
        H[self.pulse.k_ini] += 0.25* self.pulse.pulse(t) 
        H[self.pulse.k_ini-1] += 0.25* self.pulse.pulse(t)   

        return t
    
class Source:
    def __init__(self, sourcetype, t_0, s_0, k_ini):
        self.sourcetype=sourcetype
        self.t_0=t_0
        self.s_0=s_0
        self.k_ini=k_ini

    def pulse(self, time):
        
        self.time=time
        
        if self.sourcetype == 'gauss':
            pulse = np.exp(-0.5*( (self.t_0 - time) / self.s_0 )**2)
        
        return pulse