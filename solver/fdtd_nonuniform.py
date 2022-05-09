import numpy as np

class FDTD:
    def __init__(self, grid):
        self.grid = grid
        self.dualGrid = (grid[:-1] + grid[1:]) / 2.0
        speedOfLight = 1.0
        self.dx = (grid[1:] - grid[:-1])
        self.h = (self.dx[1:] + self.dx[:-1] )/2.0
        self.dt = 0.9* np.amin(self.dx)/speedOfLight 


        self.E = 0.0 * self.grid
        self.H = 0.0 * self.dualGrid

    def step(self, t):
        E = self.E
        H = self.H

        cosa = np.ones(len(E))


        E[1:-1] = E[1:-1] - (self.dt/self.h)*(H[1:] - H[:-1]) 
        H[:]    = H[:]    - (self.dt/self.dx)*(E[1:] - E[:-1]) 
        
        return t + self.dt


