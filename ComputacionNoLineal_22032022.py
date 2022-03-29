#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 13:48:47 2022

@author: juanjosefernandezmorales
"""
import numpy as np

# Funciones por separado.

def TransmissionMatrix(ep_r, sigma, d, omega, mu_r):
    ep_c  = ep_r*ep_0-1j*sigma/omega    # Permetividad electrica compleja
    mu    = mu_0*mu_r                   # Permitividad magnetica
    
    gamma = 1j*omega*np.sqrt(mu*ep_c)   # Complex propagation constant
    eta   = np.sqrt(mu/ep_c)            # Intrinsic impedance
    
    _     = gamma*d
    phi   = np.array([[np.cosh(_),           eta*np.sinh(_)], 
                      [eta**(-1)*np.sinh(_), np.cosh(_)]]) # Transmission Matrix
    return phi

def TRCoefficients(ep_r, sigma, d, omega, mu_r):
    phi = TransmissionMatrix(ep_r, sigma, d, omega, mu_r)  # Transmission Matrix
    
    _ = (phi[0,0]*eta_0 + phi[0,1] + phi[1,0]*eta_0**2 + phi[1,1]*eta_0)
    
    # Transmission coefficient
    T = 2*eta_0/_
    # Reflexion coefficient
    R = (phi[0,0]*eta_0 + phi[0,1] - phi[1,0]*eta_0**2 - phi[1,1]*eta_0)/_
    return T, R

# Creación de una clase

class Panel():
    def __init__(self, ep_r, mu_r, d, sigma, omega):
        # Datos propios del panel.
        self.ep_r  = ep_r
        self.mu_r  = mu_r
        self.d     = d
        self.sigma = sigma
        self.omega = omega
        
        self.T, self.R = TRCoefficients(self.ep_r, self.sigma, 
                                        self.d, self.omega, self.mu_r,)
    
    
    def TransmissionMatrix(ep_r, sigma, d, omega, mu_r):
        ep_c  = ep_r*ep_0-1j*sigma/omega    # Permetividad electrica compleja
        mu    = mu_0*mu_r                   # Permitividad magnetica
        
        gamma = 1j*omega*np.sqrt(mu*ep_c)   # Complex propagation constant
        eta   = np.sqrt(mu/ep_c)            # Intrinsic impedance
        
        _     = gamma*d
        phi   = np.array([[np.cosh(_),           eta*np.sinh(_)], 
                          [eta**(-1)*np.sinh(_), np.cosh(_)]]) # Transmission Matrix
        return phi

    def TRCoefficients(ep_r, sigma, d, omega, mu_r):
        phi = TransmissionMatrix(ep_r, sigma, d, omega, mu_r)  # Transmission Matrix
        
        _ = (phi[0,0]*eta_0 + phi[0,1] + phi[1,0]*eta_0**2 + phi[1,1]*eta_0)
        
        # Transmission coefficient
        T = 2*eta_0/_
        # Reflexion coefficient
        R = (phi[0,0]*eta_0 + phi[0,1] - phi[1,0]*eta_0**2 - phi[1,1]*eta_0)/_
        return T, R
    
# Datos iniciales

mu_0 = 4*np.pi*1e-7 # Permitividad magnetica del vacio
mu_r = 1            # Permitividad magnetuca relativo

ep_0 = 8.854e-12    # Permitividad electrica del vacio
ep_r = 5            # Permitividad electrica relativo

sigma = 10          # Conductividad [kS/m]
omega = 1e3*2*np.pi # Frecuencia

eta_0 = np.sqrt(mu_0/ep_0)
d     = 10e-3   # anchura


print("Matriz Transmisión:\n\n",
      TransmissionMatrix(ep_r, sigma, d, omega, mu_r), '\n')

print("Coeficientes: \n", TRCoefficients(ep_r, sigma, d, omega, mu_r), '\n')

        
my_panel = Panel(ep_r, mu_r, sigma, d, omega)
print("Coeficientes del objeto:\n", my_panel.T, my_panel.R,'\n')
