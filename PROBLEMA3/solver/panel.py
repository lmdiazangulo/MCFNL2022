#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.constants as sp

eta_0 = np.sqrt(sp.mu_0/sp.epsilon_0)

class Panel: 
    def __init__(self, thickness, epsilon_r = 1.0, sigma = 0.0, mu_r = 1.0):
        self.thickness = thickness
        self.epsilon_r = epsilon_r
        self.mu_r  = mu_r
        self.sigma = sigma

    def phi(self, omega):
        epsilon_c = self.epsilon_r*sp.epsilon_0 - complex(0,1)*self.sigma/omega
        mu_c      = self.mu_r * sp.mu_0
        gamma     = complex(0,1) * omega * np.sqrt(epsilon_c * mu_c)
            
        gd  = gamma * self.thickness
        eta = np.sqrt(mu_c / epsilon_c)
        
        phi = np.array([[np.cosh(gd),      np.sinh(gd) * eta], 
                        [np.sinh(gd) /eta, np.cosh(gd)      ]])
        
        return phi
        
    def T(self, omega):
        phi = self.phi(omega)
        den = phi[0,0]*eta_0 + phi[0,1] + phi[1,0]*eta_0**2 + phi[1,1]*eta_0
        return  2*eta_0 / den

    def R(self, omega): 
        phi = self.phi(omega)
        den = phi[0,0]*eta_0 + phi[0,1] + phi[1,0]*eta_0**2 + phi[1,1]*eta_0
        ren = phi[0,0]*eta_0 + phi[0,1] - phi[1,0]*eta_0**2 - phi[1,1]*eta_0
        return ren/den
    
def RTfunction(E1, E2, panel_min, panel_max, tEnd):
    
    E_min1 = E1[:, panel_min]
    E_max1 = E1[:, panel_max]

    E_min2 = E2[:, panel_min]
    E_max2 = E2[:, panel_max]

    R = np.abs(np.fft.fft(E_min1 - E_min2) / np.fft.fft(E_min2))
    T = np.abs(np.fft.fft(E_max1) / np.fft.fft(E_min2))
    
    freq = ((2.0*np.pi)/tEnd) * np.arange(len(E_min1))   
    
    return R, T, freq