import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
from scipy import signal
import pycwt as wavelet
import random

t = np.arange(0, 1e5, 0.1)

r = random.random()
print("r", r)

a1 = 20
a2 = 19

f1 = 1234.0
f2 = 1234.0

p1 = 0
p2 = 0

s1 = a1*np.sin(2*np.pi*f1*t + p1)
s2 = a2*np.sin(2*np.pi*f2*t + p2)

print("corr coef", np.corrcoef(s1, s2), "\n\n")

fig1 = plt.figure("signal over time")
plt.plot(t, s1, label="s1")
plt.plot(t, s2, label="s2")
plt.xlabel('Time')
plt.ylabel('Signal')
plt.grid(True)
plt.legend()


crossCorr = correlate(s1, s2, "full", "fft")
a = len(t)
b = len(crossCorr)
ratio = a/b
print("Input point count", a)
print("Output point count", b)
print("Ratio between input data length and cross-correlation is", ratio)
points = np.linspace(0, b, b)
time = points * ratio
print("Delay is", time[np.argmax(crossCorr)], "time units")

fig2 = plt.figure("corss correlation")
plt.plot(time, crossCorr)
plt.xlabel('Time')
plt.ylabel('Cross-correlation between components s1 and s2')
plt.grid(True)

fig3 = plt.figure("coherence")
f, c_xy = signal.coherence(s1, s2, nfft=1024, nperseg=1024)
f2, c_xy2 = signal.coherence(s1, s2, nfft=1024, nperseg=1024)
plt.semilogy(f, c_xy, label="scipy_signal_coherence_1")
plt.semilogy(f2, c_xy2, label="scipy_signal_coherence")
plt.xlabel('frequency')
plt.ylabel('Coherence')

plt.legend()
plt.grid(True)
plt.show()


fs = 10e3
print("fs", fs)
N = 1e5
amp = 20
freq = 1234.0
noise_power = 0.001 * fs / 2
time = np.arange(N) / fs
b, a = signal.butter(2, 0.25, 'low')
x = np.random.normal(scale=np.sqrt(noise_power), size=time.shape)
y = signal.lfilter(b, a, x)
x += amp*np.sin(2*np.pi*freq*time)
y += np.random.normal(scale=0.1*np.sqrt(noise_power), size=time.shape)

f, Cxy = signal.coherence(x, y, fs, nperseg=1024)
fig4 = plt.figure("coh 2")
plt.semilogy(f, Cxy)
plt.xlabel('frequency [Hz]')
plt.ylabel('Coherence')

fig5 = plt.figure("coh time")
plt.plot(x)
plt.plot(y)
plt.xlabel('time')
plt.show()




