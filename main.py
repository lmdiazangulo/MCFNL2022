import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as npf

tIni = 0
tEnd = 500e-6
N = int(1e5 + 1)

s0 = 1e-6
t0 = 10*s0

t = np.linspace(tIni, tEnd, N)
g = np.exp(- np.power(t - t0, 2)/(2 * s0**2))

fq = npf.fftshift(npf.fftfreq(N, t[1]-t[0]))

plt.plot(fq, npf.fftshift(np.abs(npf.fft(g))))
plt.xlim(-100e3, 100e3)

plt.show()

print("END")