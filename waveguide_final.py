import numpy as np
import scipy.io

# Grid x, y
N = 15000 #Número de simulaciones en el tiempo
nx, ny = 101, 101 #Número de puntos en el grid
xs, dx = np.linspace(0, 1, nx+1, endpoint=True, retstep=True)
ys, dy = np.linspace(0, 3, ny+1, endpoint=True, retstep=True)


x, y = np.meshgrid(0.5*(xs[:-1]+xs[1:]), 0.5*(ys[:-1]+ys[1:]))

# Inicializa E_y, E_x en el primer tiempo
e_x = np.zeros((N, nx, ny+1), dtype=np.double)  
e_y = np.zeros((N, nx+1, ny), dtype=np.double) 

c = 1

dt = 0.01/(c*np.sqrt(1/dx**2+1/dy**2))


# Inicializa H_z en el primer tiempo
h_z = np.zeros((N, nx, ny), dtype=np.double)  


#hz_0 = np.exp(-(y-0.1)**2/0.035**2-(x-0.5)**2/0.5**2)
#hz_0 = np.exp(-(y-0.1)**2/0.05**2+0*x)
hz_0 = np.exp(-(y-0.1)**2/0.1**2-(x-0.35)**2/0.1**2)+np.exp(-(y-0.1)**2/0.1**2-(x-0.65)**2/0.1**2)

h_z[0,:,:] = hz_0

# Updates
for t in range(N-1):
    
    e_x[t+1, :, 1:-1] = e_x[t, :, 1:-1] + dt/dy*(h_z[t, :, 1:]-h_z[t, :, :-1])
    e_y[t+1, 1:-1, :] = e_y[t, 1:-1, :] - dt/dx*(h_z[t, 1:, :]-h_z[t, :-1, :]) 

    e_x[t+1, :, 0] = 0 #Pared de izquierda (ejes)
    e_x[t+1, :, -1] = 0 #Pared de derecha 
    
    #e_y[t+1, -1, :] = 0 #Pared de arriba
    #e_y[t+1, 0, :] = 0 #Pared de abajo
    e_y[t+1, 0, :] = e_y[t, 1, :]+(dt-dx)/(dt+dx)*(e_y[t+1, 1, :]-e_y[t, 0, :])
    e_y[t+1, -1, :] = e_y[t, -2, :]+(dt-dx)/(dt+dx)*(e_y[t+1, -2, :]-e_y[t, -1, :]) 


    h_z[t+1, :, :] = h_z[t, :, :] + dt/dy*(e_x[t, :, 1:]-e_x[t, :, :-1]) - dt/dx*(e_y[t, 1:, :]-e_y[t, :-1, :])


# ------------------------------------------------------------------------------------------

scipy.io.savemat('h_z_temp.mat', mdict={'out': h_z}, oned_as='row')
