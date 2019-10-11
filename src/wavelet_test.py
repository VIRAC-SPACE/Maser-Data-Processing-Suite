#! /usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from obspy.imaging.cm import obspy_sequential
from obspy.signal.tf_misfit import cwt
import argparse
import pycwt as wavelet


from help import *
from parsers._configparser import ConfigParser


def parseArguments():
    parser = argparse.ArgumentParser(description='''Wavelet. ''', epilog="""Wavelet.""")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 2.0')
    args = parser.parse_args()
    return args


def get_args(key):
    args = parseArguments()
    return str(args.__dict__[key])


def get_configs(key, value):
    config_file_path = get_args("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(config_file_path)
    return config.getConfig(key, value)


source = "cepa"
component_count = len(get_configs("velocities", source).replace(" ", "").split(","))
file = "/home/janis/Documents/maser/monitoring/cepa_6668.txt"
data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file), component_count + 1))
a = 5
x = correctNumpyReadData(data[:, [a]])
time = correctNumpyReadData(data[:, [0]])

t0 = time[0]
tmax = time[-1]
dt = (tmax - t0)/len(time)
dt = dt / (24 * 60 * 60)

scalogram = cwt(x, dt , 8, np.min(x), np.max(x), wl="morlet")

fig = plt.figure()
ax = fig.add_subplot(111)

x, y = np.meshgrid(time, np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), scalogram.shape[0]))

ax.pcolormesh(x, y, np.abs(scalogram), cmap=obspy_sequential)
ax.set_xlabel("Time")
ax.set_ylabel("Frequency [Hz]")
ax.set_yscale('log')
#ax.set_ylim(np.min(x), np.max(x))
plt.show()

x = correctNumpyReadData(data[:, [a]])
dt = (tmax - t0)/len(time)
N = len(x)
p = np.polyfit(time - t0, x, 1)
dat_notrend = x - np.polyval(p, time - t0)
std = dat_notrend.std()
var = std ** 2
dat_norm = dat_notrend / std

mother = wavelet.Morlet(6)
s0 = 2 * dt
dj = 1 / 12
J = 7 / dj
alpha, _, _ = wavelet.ar1(x)  # Lag-1 autocorrelation for red noise

wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dat_norm, dt, dj, s0, J, mother)
iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std

power = (np.abs(wave)) ** 2
fft_power = np.abs(fft) ** 2
period = 1 / freqs
power /= scales[:, None]

signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha, significance_level=0.9999, wavelet=mother)
sig95 = np.ones([1, N]) * signif[:, None]
sig95 = power / sig95

glbl_power = power.mean(axis=1)
dof = N - scales
glbl_signif, tmp = wavelet.significance(var, dt, scales, 1, alpha, significance_level=0.95, dof=dof, wavelet=mother)

'''
sel = find((period >= 2) & (period < 8))
Cdelta = mother.cdelta
scale_avg = (scales * np.ones((N, 1))).transpose()
scale_avg = power / scale_avg  # As in Torrence and Compo (1998) equation 24
scale_avg = var * dj * dt / Cdelta * scale_avg[sel, :].sum(axis=0)
scale_avg_signif, tmp = wavelet.significance(var, dt, scales, 2, alpha, significance_level=0.95, dof=[scales[sel[0]], scales[sel[-1]]], wavelet=mother)
'''

# Prepare the figure
plt.close('all')
plt.ioff()
figprops = dict(figsize=(11, 8), dpi=72)
fig = plt.figure(**figprops)

# First sub-plot, the original time series anomaly and inverse wavelet
# transform.
ax = plt.axes([0.1, 0.75, 0.65, 0.2])
ax.plot(time, iwave, '-', linewidth=1, color=[0.5, 0.5, 0.5], label="iwave")
ax.plot(time, x, 'k', linewidth=1.5, label="orginal data")
ax.legend()

# Second sub-plot, the normalized wavelet power spectrum and significance
# level contour lines and cone of influece hatched area. Note that period
# scale is logarithmic.
bx = plt.axes([0.1, 0.37, 0.65, 0.28], sharex=ax)
levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
bx.contourf(time, np.log2(period), np.log2(power), np.log2(levels), extend='both', cmap=plt.cm.viridis)
extent = [np.min(time), np.max(time), 0, max(period)]
bx.contour(time, np.log2(period), sig95, [-99, 1], colors='k', linewidths=2, extent=extent)
bx.fill(np.concatenate([time, time[-1:] + dt, time[-1:] + dt, time[:1] - dt, time[:1] - dt]), np.concatenate([np.log2(coi), [1e-9], np.log2(period[-1:]), np.log2(period[-1:]), [1e-9]]), 'k', alpha=0.3, hatch='x')
#bx.set_title('b) {} Wavelet Power Spectrum ({})'.format(label, mother.name))
bx.set_ylabel('Period (years)')
#
Yticks = 2 ** np.arange(np.ceil(np.log2(period.min())), np.ceil(np.log2(period.max())))
bx.set_yticks(np.log2(Yticks))
bx.set_yticklabels(Yticks)

# Third sub-plot, the global wavelet and Fourier power spectra and theoretical
# noise spectra. Note that period scale is logarithmic.
cx = plt.axes([0.77, 0.37, 0.2, 0.28], sharey=bx)
cx.plot(glbl_signif, np.log2(period), 'k--')
cx.plot(var * fft_theor, np.log2(period), '--', color='#cccccc')
cx.plot(var * fft_power, np.log2(1./fftfreqs), '-', color='#cccccc', linewidth=1.)
cx.plot(var * glbl_power, np.log2(period), 'k-', linewidth=1.5)
cx.set_title('c) Global Wavelet Spectrum')
#cx.set_xlabel(r'Power [({})^2]'.format(units))
cx.set_xlim([0, glbl_power.max() + var])
cx.set_ylim(np.log2([period.min(), period.max()]))
cx.set_yticks(np.log2(Yticks))
cx.set_yticklabels(Yticks)
plt.setp(cx.get_yticklabels(), visible=False)


plt.show()