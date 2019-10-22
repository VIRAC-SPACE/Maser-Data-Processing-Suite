import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
from scipy import signal

fs = 1000
t = []

ti = 0
while ti < 1:
    t += [ti]
    ti += 1/fs
t = np.array(t)

#print("time", t)

a1 = 15000
a2 = 15000
a3 = 15000

f1 = 1234.0
f2 = f1 * 1.5
f3 = f2 * 2

p1 = 0
p2 = 0
p3 = 0

s1 = a1*np.sin(2*np.pi*f1*t + p1)
s2 = a2*np.sin(2*np.pi*f2*t + p2)
s3 = a3*np.sin(2*np.pi*f3*t + p3)

#print("corr coef", np.corrcoef(s1, s2), "\n\n")

fig1 = plt.figure("signal over time")
plt.plot(t, s1, label="s1")
plt.plot(t, s2, label="s2")
plt.plot(t, s3, label="s2")
plt.xlabel('Time')
plt.ylabel('Signal')
plt.grid(True)
plt.legend()
plt.show()

signal_ft_s1 = np.fft.fft(s1)
signal_ft_s2 = np.fft.fft(s2)
signal_ft_s3 = np.fft.fft(s3)

dt_s1 = 1/f1
dt_s2 = 1/f2
dt_s3 = 1/f3

ftfreqs_s1 = 2 * np.pi * np.fft.fftfreq(len(signal_ft_s1), dt_s1)
ftfreqs_s2 = 2 * np.pi * np.fft.fftfreq(len(signal_ft_s2), dt_s2)
ftfreqs_s3 = 2 * np.pi * np.fft.fftfreq(len(signal_ft_s3), dt_s2)

f0 = 6
s0 = 2
dj = 1 / 12
J = -1
n0 = len(s1)
J = np.int(np.round(np.log2(n0 * dt_s1/ s0) / dj))
J = abs(J)
sj = s0 * 2 ** (np.arange(0, J + 1) * dj)
sj_col = sj[:, np.newaxis]

f_s1 = sj_col * ftfreqs_s1
f_s2 = sj_col * ftfreqs_s2
f_s3 = sj_col * ftfreqs_s3

wavelet_psi_ft_s1 = (np.pi ** -0.25) * np.exp(-0.5 * (f_s1 - f0) ** 2)
wavelet_psi_ft_s2 = (np.pi ** -0.25) * np.exp(-0.5 * (f_s2 - f0) ** 2)
wavelet_psi_ft_s3 = (np.pi ** -0.25) * np.exp(-0.5 * (f_s3 - f0) ** 2)

psi_ft_bar_s1 = ((sj_col * ftfreqs_s1[1] * len(signal_ft_s1)) ** .5 * np.conjugate(wavelet_psi_ft_s1 ))
psi_ft_bar_s2 = ((sj_col * ftfreqs_s2[1] * len(signal_ft_s2)) ** .5 * np.conjugate(wavelet_psi_ft_s2 ))
psi_ft_bar_s3 = ((sj_col * ftfreqs_s3[1] * len(signal_ft_s3)) ** .5 * np.conjugate(wavelet_psi_ft_s3 ))

w_s1 = np.fft.ifft(signal_ft_s1 * psi_ft_bar_s1, axis=1)
w_s2 = np.fft.ifft(signal_ft_s2 * psi_ft_bar_s2, axis=1)
w_s3 = np.fft.ifft(signal_ft_s3 * psi_ft_bar_s3, axis=1)

sel_s1 = np.invert(np.isnan(w_s1).all(axis=1))
sel_s2 = np.invert(np.isnan(w_s2).all(axis=1))
sel_s3 = np.invert(np.isnan(w_s3).all(axis=1))

wavelet_flambda = (4 * np.pi) / (f0 + np.sqrt(2 + f0 ** 2))

freqs_s1 = 1 / (wavelet_flambda * sj)
freqs_s2 = 1 / (wavelet_flambda * sj)
freqs_s3 = 1 / (wavelet_flambda * sj)

if np.any(sel_s1):
        sj = sj[sel_s1]
        freqs_s1 = freqs_s1[sel_s1]
        w_s1 = w_s1[sel_s1, :]

if np.any(sel_s2):
        sj = sj[sel_s2]
        freqs_s2 = freqs_s3[sel_s2]
        w_s2 = w_s2[sel_s2, :]

if np.any(sel_s3):
        sj = sj[sel_s3]
        freqs_s3 = freqs_s3[sel_s3]
        w_s3 = w_s3[sel_s3, :]

w_s1 = w_s1[:, :n0]
w_s2 = w_s2[:, :n0]
w_s3 = w_s3[:, :n0]

#print(w_s1, w_s2, w_s3)

power_s1 = (np.abs(w_s1)) ** 2
power_s2 = (np.abs(w_s2)) ** 2
power_s3 = (np.abs(w_s3)) ** 2

period_s1 = 1/freqs_s1
period_s2 = 1/freqs_s2
period_s3 = 1/freqs_s2

print(power_s1)

levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]

plt.subplot(231)
plt.xlabel('Time')
plt.ylabel('Period s1')
plt.contourf(t, np.log2(period_s1), np.log2(power_s1), np.log2(levels), extend='both', cmap=plt.cm.viridis)

plt.subplot(232)
plt.xlabel('Time')
plt.ylabel('Period s2')
plt.contourf(t, np.log2(period_s2), np.log2(power_s2), np.log2(levels), extend='both', cmap=plt.cm.viridis)

plt.subplot(233)
plt.xlabel('Time')
plt.ylabel('Period s2')
plt.contourf(t, np.log2(period_s3), np.log2(power_s3), np.log2(levels), extend='both', cmap=plt.cm.viridis)

# bez log
plt.subplot(234)
plt.xlabel('Time')
plt.ylabel('Period s1')
plt.contourf(t, period_s1, power_s1, levels, extend='both', cmap=plt.cm.viridis)

plt.subplot(235)
plt.xlabel('Time')
plt.ylabel('Period s2')
plt.contourf(t, period_s2, power_s2, levels, extend='both', cmap=plt.cm.viridis)

plt.subplot(236)
plt.xlabel('Time')
plt.ylabel('Period s2')
plt.contourf(t, period_s3, power_s3, levels, extend='both', cmap=plt.cm.viridis)

plt.show()


'''
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
f, c_xy = signal.coherence(s1, s2, fs=fs)
plt.plot(f, c_xy, label="s1 and s2")

f, c_xy = signal.coherence(s1, s3, fs=fs)
plt.plot(f, c_xy, label="s1 and s3")

f, c_xy = signal.coherence(s3, s2, fs=fs)
plt.plot(f, c_xy, label="s3 and s2")

plt.xlabel('frequency')
plt.ylabel('Coherence')

plt.legend()
plt.grid(True)
plt.show()

'''



