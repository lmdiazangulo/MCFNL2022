#!/usr/bin/env python3

from material.panel import *

epsilonr = 5
sigma = 10
d = 10e-03
w = 2*np.pi*1000#entre 1khz y 1GHz

print("Coeficiente de reflexión:", Panel(epsilonr,sigma,d,w).Reflexion())
print("Coeficiente de transmisión:", Panel(epsilonr,sigma,d,w).Transmision())


        
        

    
