import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
from scipy import signal
import pycwt as wavelet

t = np.arange(0, 100, 0.1)

a1 = 1
a2 = 2

f1 = 0.5
f2 = 0.6

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
f, c_xy = signal.coherence(s1, s2, nfft=256)
plt.semilogy(f, c_xy, label="scipy_signal_coherence")
plt.xlabel('frequency')
plt.ylabel('Coherence')
plt.cohere(s1, s2, NFFT=256, label="pyplot_cohere")
plt.legend()
plt.grid(True)
plt.show()




