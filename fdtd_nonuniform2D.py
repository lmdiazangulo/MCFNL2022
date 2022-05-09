import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation


xMax = 1
xMin = 0

#Create non-uniform grid
#x = xMin
#v = [x]
#step = 0.005
#while x < xMax:
    #x = x+step
    #v.append(x)
    #step = step/0.97
    
#vx = np.array(v)
#vy = np.copy(vx)

a = np.linspace(xMin, xMax/3, 20,  endpoint = False)
b = np.linspace(xMax/3, 2*xMax/3, 40, endpoint = False)
c = np.linspace(2*xMax/3, xMax, 20)
d = np.concatenate((a,b,c))

vx = d
vy = np.copy(vx)


#vx = np.array(v)
#vy = np.copy(vx)


#vx = np.linspace(xMin, xMax, 100)
#vy = np.copy(vx)


dx = (vx[1:] - vx[:-1])
dy = (vy[1:] - vy[:-1])
hx = (dx[1:] + dx[:-1])/2.0
hy = (dy[1:] + dy[:-1])/2.0

#Initialize Ex, Ey 
Ex = np.zeros((len(vx)-1, len(vy)), dtype=np.double)  
Ey = np.zeros((len(vx), len(vy)-1), dtype=np.double) 

speedOfLight = 1.0
dt = 0.2*(1/np.sqrt((1/np.amin(dx)**2 + 1/np.amin(dy)**2)))/speedOfLight     # time step

# Create grid 
x, y = np.meshgrid(0.5*(vx[:-1]+vx[1:]), 0.5*(vy[:-1]+vy[1:]))

fig = plt.figure()
ax = axes3d.Axes3D(fig)
plt.rcParams['figure.figsize'] = (9,7)

mean = (xMax + xMin)/2
spread = 0.00005

Hz = np.exp( - (np.power(x - mean,2)/2.0/np.sqrt((spread)) + np.power(y - mean,2)/2.0/np.sqrt((spread))) )  

def faraday(Ex, Ey, Hz) : 
    "faraday equation Bz(t+dt/2) -> Bz(t-dt/2) + dt f(E(t))"
    return Hz + ((dt / dy)*(Ex[:, 1:]-Ex[:, :-1]) - (dt / dx)*(Ey[1:, :]-Ey[:-1, :]))

    
def ampere_maxwell(Hz, Ex, Ey):
    " Ampere-Maxwell equation E(t+dt) -> E(t) + dt g(Bz(t+dt/2)) "
    Ex[:, 1:-1] += (dt / hy)*(Hz[:, 1:]-Hz[:, :-1])
    Ey[1:-1, :] += - (dt / hx[...,None])*(Hz[1:, :]-Hz[:-1, :]) 
    return Ex, Ey


def update(i, ax, fig):
    ax.cla()

    global Ex, Ey, Hz

    for j in range(10):
        Hz = faraday(Ex, Ey, Hz)
        Ex, Ey = ampere_maxwell(Hz, Ex, Ey)
    
    wframe = ax.plot_wireframe(x, y, Hz, rstride=2, cstride=2, color='indigo')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('$E_z$')
    ax.set_zlim(-1, 1)
    return wframe,


ani = animation.FuncAnimation(fig, update,
                              frames=range(100),
                              fargs=(ax, fig), interval=200, blit=True)
