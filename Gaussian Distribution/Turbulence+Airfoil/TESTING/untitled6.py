from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

T = 10  # Duration in seconds
f0 = 100  # Fundamental frequency
Fs = 1000  # Sampling frequency

# Time domain signal
t = np.arange(0, T*Fs)/Fs
print len(t)
x = np.sin(2*np.pi*f0*t)
N = x.size

# DFT
X = np.fft.fft(x)
X_db = 20*np.log10(2*np.abs(X)/N)
#f = np.fft.fftfreq(N, 1/Fs)
f = np.arange(0, N)*Fs/N

plt.plot(f, X_db)
#plt.xlim(0,1000)
plt.grid()
plt.show()