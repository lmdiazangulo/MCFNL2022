import numpy as np
import math
import scipy.constants
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ==== Preamble ===============================================================
c0   = scipy.constants.speed_of_light
mu0  = scipy.constants.mu_0
eps0 = scipy.constants.epsilon_0
imp0 = math.sqrt(mu0 / eps0)
pi = math.pi

# ==== Inputs / Pre-processing ================================================ 
# ---- Problem definition -----------------------------------------------------
L         = 10.0
dx        = 0.1
dy        = 0.1
finalTime = 0.5*10**(-7)
cfl       = .99

gridEX = np.linspace(0,      L,        num=int(L/dx + 1), endpoint=True)
gridEY = np.linspace(0,      L,        num=int(L/dy + 1), endpoint=True)
gridHX = np.linspace(dx/2.0, L-dx/2.0, num=int(L/dx),   endpoint=True)
gridHY = np.linspace(dy/2.0, L-dy/2.0, num=int(L/dy),   endpoint=True)

# ---- Materials --------------------------------------------------------------
eps = np.zeros((gridEX.size, gridEY.size), dtype = complex)

epsr = np.zeros((gridEX.size, gridEY.size))
sigma = np.zeros((gridEX.size, gridEY.size))

for i1 in range(10,91):
    for j1 in range(10,91):
        epsr[i1][j1] = 1.0
        sigma[i1][j1] = 0.0
        
for i2 in range(0,10):
    epsr[i2][:] = 2.0
    sigma[i2][:] = 5.0
    
for ix in range(gridEX.size):  
    for j2 in range(0,10):
        epsr[ix][j2] = 2.0
        sigma[ix][j2] = 5.0
        
for i3 in range(91,gridEX.size): 
    for jy in range(gridEX.size):
        epsr[jy][i3] = 2.0
        sigma[jy][i3] = 5.0
    
for j3 in range(91,gridEY.size):
    epsr[:][j3] = 2.0
    sigma[:][j3] = 5.0

eps[:][:] = epsr[:][:] + complex(0,-1)*sigma[:][:]

# ---- Sources ----------------------------------------------------------------
# Initial field H

spread   = 0.1
ampl = 0.1

# Frequecy Range in Physical Anechoic Chambers : 200 Hz - 20 kHz
f = 10**(3)

center   = (L/2.0, L/2.0)

initialH = np.zeros((gridHX.size, gridHY.size))
for i in range(gridHX.size):
    for j in range(gridHY.size):
        initialH[i,j] = math.exp( 
            - ((gridHX[i]-center[0])**2 + (gridHY[j]-center[1])**2) /
            math.sqrt(2.0) / spread)
        
        # initialH[i,j] = ampl*np.cos((2*pi*f)*((gridHX[i]-center[0]) +  \
                                              # (gridHY[j]-center[1]))) 
        
 
# ---- Output requests --------------------------------------------------------
samplingPeriod = 0.0
 
# ==== Processing =============================================================
# ---- Solver initialization --------------------------------------------------
dt = cfl * dx / c0 / np.sqrt(2)
numberOfTimeSteps = int( finalTime / dt )

if samplingPeriod == 0.0:
    samplingPeriod = dt 
nSamples  = int( math.floor(finalTime/samplingPeriod) )

probeH    = np.zeros((gridHX.size, gridHY.size, nSamples))
probeTime = np.zeros(nSamples) 

# probeE    = np.zeros((gridEX.size, gridEY.size, nSamples))
# probeTime = np.zeros(nSamples) 

exOld = np.zeros((gridEX.size, gridEY.size))
exNew = np.zeros((gridEX.size, gridEY.size))
eyOld = np.zeros((gridEX.size, gridEY.size))
eyNew = np.zeros((gridEX.size, gridEY.size))
hzOld = np.zeros((gridHX.size, gridHY.size))
hzNew = np.zeros((gridHX.size, gridHY.size))

if 'initialH' in locals():
    hzOld = initialH

# Determines recursion coefficients
cEx = np.zeros((gridEX.size, gridEY.size), dtype = complex) 
cEy = np.zeros((gridEX.size, gridEY.size), dtype = complex)

for i in range(gridEX.size):
    for j in range(gridEY.size):
        cEx[i][j] = dt / (eps0*(eps[i][j] + sigma[i][j]*dt/2)) / dx
        cEy[i][j] = dt / (eps0*(eps[i][j] + sigma[i][j]*dt/2)) / dy

cHx = dt / mu0 / dx
cHy = dt / mu0 / dy

# ---- Time integration -------------------------------------------------------
print('--- FDTD 2D ---')
tic = time.time();

t = 0.0
for n in range(numberOfTimeSteps):
    # --- Updates E field ---
    for i in range(1, gridEX.size-1):
        for j in range(1, gridEY.size-1):
            exNew[i][j] = ((eps[i][j] - sigma[i][j]*dt/2)/(eps[i][j] + sigma[i][j]*dt/2))*exOld[i][j] \
                + cEy[i][j] * (hzOld[i][j] - hzOld[i  ][j-1])
            eyNew[i][j] = ((eps[i][j] - sigma[i][j]*dt/2)/(eps[i][j] + sigma[i][j]*dt/2))*eyOld[i][j] \
                - cEx[i][j] * (hzOld[i][j] - hzOld[i-1][j  ])
     
    # PEC
    exNew[ :][ 0] = 0.0
    exNew[ :][-1] = 0.0
    eyNew[ 0][ :] = 0.0
    eyNew[-1][ :] = 0.0 
    

    # --- Updates H field ---
    for i in range(gridHX.size):
        for j in range(gridHX.size):
            hzNew[i][j] = hzOld[i][j] - cHx * (eyNew[i+1][j  ] - eyNew[i][j]) +\
                                        cHy * (exNew[i  ][j+1] - exNew[i][j])
    
    # PMC
    hzNew[ :][ 0] = 0.0
    hzNew[ :][-1] = 0.0
    for i in range(gridHX.size):
        hzNew[ i][0] = 0.0
        hzNew[ i][-1] = 0.0
    
        
   
    # --- Updates output requests ---
    probeH[:,:,n] = hzNew[:,:]
    probeTime[n] = t
    
    # --- Updates fields and time 
    exOld[:] = exNew[:]
    eyOld[:] = eyNew[:]
    hzOld[:] = hzNew[:]
    t += dt
    # print ("Time step: %d of %d" % (n, numberOfTimeSteps-1))

tictoc = time.time() - tic;
print('--- End of simulation ---')
print("CPU Time: %f [s]" % tictoc)

# ==== Post-processing ========================================================

        
# --- Creates animation ---
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(1, 1, 1)
#ax = plt.axes(xlim=(gridE[0], gridE[-1]), ylim=(-1.1, 1.1))
ax.set_xlabel('X coordinate ')
ax.set_ylabel('Y coordinate ')
line = plt.imshow(probeH[:,:,0], animated=True, vmin=-0.15, vmax=0.15, cmap ="RdBu")
timeText = ax.text(0.02, 0.95, '', transform=ax.transAxes)
rectangle1 = plt.Rectangle((9 ,8), 81, 82, edgecolor="black", fill=False)
plt.gca().add_patch(rectangle1)

def init():
    line.set_array(probeH[:,:,0])
    timeText.set_text('')
    return line, timeText

def animate(i):
    line.set_array(probeH[:,:,i])
    timeText.set_text('Time = %2.1f [ns]' % (probeTime[i]*1e9))
    return line, timeText

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=nSamples, interval=50, blit=True)

anim.save('simulacion2d.gif', fps = 30)
plt.show()

print('=== END ===')