import numpy as np

class Panel:
    """
    Define un panel según su permitividad relativa, su conductancia y grosor
    """
    
    def __init__(self,epsilonr,sigma,d,w):
        self.epsilonr = epsilonr
        self.sigma = sigma
        self.d = d
        self.w = w
    
    def TransmissionMatrix(self):
        epsilon0 = 8.8542e-12
        mu = 4*np.pi-7
        epsilonc = self.epsilonr*epsilon0 - 1j*self.sigma/self.w
        gamma = 1j*self.w*np.sqrt(mu*epsilonc)
        eta = np.sqrt(mu/epsilonc)
    
        #matriz de transmisión
        phi = np.array([[np.cosh(gamma*self.d),eta*np.sinh(gamma*self.d)],[(1/eta)*np.sinh(gamma*self.d),np.cosh(gamma*self.d)]])
    
        return phi
      
    def Reflexion(self):
        phi = self.TransmissionMatrix()
        eta0 = 376.73031346177
        R = (phi[0][0]*eta0 + phi[0][1] - phi[1][0]*eta0**2 - phi[1][1]*eta0)/(phi[0][0]*eta0 + phi[0][1] + phi[1][0]*eta0**2 + phi[1][1]*eta0)
        return R
    
    def Transmision(self):
        phi = self.TransmissionMatrix()
        eta0 = 376.73031346177
        T = 2*eta0/(phi[0][0]*eta0 + phi[0][1] + phi[1][0]*eta0**2 + phi[1][1]*eta0)
        return T