#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
from scipy.signal import correlate, resample
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
import pycwt as wavelet

from help import *
from parsers._configparser import ConfigParser

def parseArguments():
    parser = argparse.ArgumentParser(description='''Cross - Correlatet two maser componets. ''', epilog="""CrossCorr.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("file", help="Experiment source", type=str, default="")
    parser.add_argument("a", help="componet A ", type=int)
    parser.add_argument("b", help="componet B", type=str)
    parser.add_argument("--index1", help="index 1", type=str, default="")
    parser.add_argument("--index2", help="index 2", type=str, default="")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 2.0')
    args = parser.parse_args()
    return args

def getArgs(key):
    return str(parseArguments().__dict__[key])

def getConfigs(key, value):
    configFilePath = getArgs("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)

def getData(file):
    a = int(getArgs("a"))
    b = int(getArgs("b"))
    compunetCount = len(getConfigs("velocities", getArgs("source")).replace(" ", "").split(","))
    data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file),compunetCount + 1))
    x = correctNumpyReadData(data[:, [a]])
    y = correctNumpyReadData(data[:, [b]])
    time = correctNumpyReadData(data[:, [0]])
    return (x, y, time)

def main():
    index1 = getArgs("index1")
    index2 = getArgs("index2")
    file = getConfigs("paths", "monitoringFilePath") + "/" + getArgs("file")

    if index1 == "":
        index1 = 0

    if index2 == "":
        index2 = len(getData(file)[0])

    index1 = int(index1)
    index2 = int(index2)

    x = getData(file)[0][index1:index2]
    y = getData(file)[1][index1:index2]
    time = getData(file)[2][index1:index2]

    print("Resampling\n")

    xResample = resample(x, len(x), t=time)
    yResample = resample(y, len(y), t=time)
    print("len of time", len(time), len(xResample[0]))

    plt.subplot(1, 2, 1)
    plt.plot(xResample[1], xResample[0])
    plt.xlabel('Time')
    plt.ylabel('Resampled components A' )
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(yResample[1], yResample[0])
    plt.xlabel('Time')
    plt.ylabel('Resampled components B')
    plt.grid(True)

    plt.show()

    print("corr coef befor resampling", np.corrcoef(x, y), "\n\n")
    print ("corr coef after resampling", np.corrcoef(xResample[0], yResample[0]), "\n\n")

    print("Direct")
    crossCorr = np.correlate(xResample[0], yResample[0], "full")
    a = len(xResample[0])
    b = len(crossCorr)
    ratio = a/b
    print("Input point count", a)
    print("Output point count", b)
    print("Ratio between input data length and cross-correlation is", ratio)
    points = np.linspace(0, b, b)
    time = points * ratio
    print("Delay is", time[np.argmax(crossCorr)], "time units")

    plt.subplot(1,2,1)
    plt.plot(time, crossCorr)
    plt.xlabel('Time')
    plt.ylabel('Cross-correlation between components ' + getArgs("a") + " " + getArgs("b"))
    plt.grid(True)
    plt.title("Direct")

    print("\n\nFFT")
    crossCorr = correlate(xResample[0], yResample[0], "full", "fft")
    a = len(x)
    b = len(crossCorr)
    ratio = a / b
    print("Input point count", a)
    print("Output point count", b)
    print("Ratio between input data length and cross-correlation is", ratio)
    points = np.linspace(0, b, b)
    time = points * ratio
    print("Delay is", time[np.argmax(crossCorr)], "time units")

    plt.subplot(1, 2, 2)
    plt.plot(time, crossCorr)
    plt.xlabel('Time')
    plt.ylabel('Cross-correlation between components ' + getArgs("a") + " " + getArgs("b"))
    plt.grid(True)
    plt.title("FFT")
    plt.show()

    time = xResample[1]  #getData(file)[2][index1:index2]
    dt = np.diff(time)[0]
    fs = 1/dt
    print("dt, fs", dt, fs)

    f, c_xy = signal.coherence(xResample[0], yResample[0], fs=fs/(24*60*60), nfft=len(time))
    plt.plot(f, c_xy, label="coherence")
    plt.xlabel('frequency')
    plt.ylabel('Coherence')
    plt.legend()
    plt.grid(True)
    plt.show()

    n = len(time)
    mother = wavelet.Morlet(6)
    slevel = 0.95  # Significance level
    dj = 1 / 12  # Twelve sub-octaves per octaves
    s0 = 2  # 2 * dt                   # Starting scale, here 6 months
    J = -1  # 7 / dj

    p_a, residuals_a, rank_a, singular_values_a, rcond_a = np.polyfit(time, xResample[0], 1, full=True)
    print("residuals_a, rank_a, singular_values_a, rcond_a",residuals_a, rank_a, singular_values_a, rcond_a)
    dat_notrend_a = xResample[0] - np.polyval(p_a, time)
    std_a = dat_notrend_a.std()
    dat_norm_a = dat_notrend_a / std_a

    p_b, residuals_b, rank_b, singular_values_b, rcond_b = np.polyfit(time, yResample[0], 1, full=True)
    print("residuals_b, rank_b, singular_values_b, rcond_b", residuals_b, rank_b, singular_values_b, rcond_b)
    dat_notrend_b = yResample[0] - np.polyval(p_b, time)
    std_b = dat_notrend_b.std()
    dat_norm_b = dat_notrend_b / std_b

    fig = plt.figure("Normal data")
    plt.plot(time, dat_norm_a, label="dat_norm_a")
    plt.plot(time, dat_norm_b, label="dat_norm_b")
    plt.xlabel('Time')
    plt.ylabel('Flux')
    plt.legend()
    plt.grid(True)
    plt.show()

    wave_a, scales_a, freqs_a, coi_a, fft_a, fftfreqs_a = wavelet.cwt(dat_norm_a, dt, dj, s0, J, mother)
    signif_a, fft_theor_a = wavelet.significance(1.0, dt, scales_a, 0, 0.0, significance_level=slevel, wavelet=mother)
    power_a = (np.abs(wave_a)) ** 2
    fft_power_a = np.abs(fft_a) ** 2
    period_a = 1 / freqs_a

    wave_b, scales_b, freqs_b, coi_b, fft_b, fftfreqs_b = wavelet.cwt(dat_norm_b, dt, dj, s0, J, mother)
    signif_b, fft_theor_b = wavelet.significance(1.0, dt, scales_b, 0, 0.0, significance_level=slevel, wavelet=mother)
    power_b = (np.abs(wave_b)) ** 2
    fft_power_b = np.abs(fft_b) ** 2
    period_b = 1 / freqs_b

    levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]

    plt.figure("Continuous wavelet power spectrum")

    plt.subplot(121)
    plt.xlabel('Time')
    plt.ylabel('Period for components ' + getArgs("a"))
    plt.contourf(time, np.log2(period_a), np.log2(power_a), np.log2(levels), extend='both', cmap=plt.cm.viridis)

    plt.subplot(122)
    plt.xlabel('Time')
    plt.ylabel('Period for components ' + getArgs("b"))
    plt.contourf(time, np.log2(period_b), np.log2(power_b), np.log2(levels), extend='both', cmap=plt.cm.viridis)
    plt.show()

    W12, cross_coi, freq, signif = wavelet.xwt(dat_norm_a, dat_norm_b, dt, dj=dj, s0=-1, J=-1, significance_level=0.8646, wavelet='morlet', normalize=True)

    cross_power = np.abs(W12) ** 2
    cross_sig = np.ones([1, n]) * signif[:, None]
    cross_sig = cross_power / cross_sig  # Power is significant where ratio > 1
    cross_period = 1 / freq

    plt.figure("Cross wavelet transform power")
    plt.xlabel('Time')
    plt.ylabel('Power')
    plt.contourf(time, np.log2(cross_period), np.log2(cross_power), np.log2(levels), extend='both', cmap=plt.cm.viridis)
    plt.show()

    WCT, aWCT, corr_coi, freq, sig = wavelet.wct(dat_norm_a, dat_norm_b, dt, dj=dj, s0=-1, J=-1, significance_level=0.8646, wavelet='morlet', normalize=True, cache=False)

    cor_sig = np.ones([1, n]) * sig[:, None]
    cor_sig = np.abs(WCT) / cor_sig  # Power is significant where ratio > 1
    cor_period = 1 / freq

    angle = 0.5 * np.pi - aWCT
    u, v = np.cos(angle), np.sin(angle)

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.55, 0.05, 0.35])
    cbar_ax_1 = fig.add_axes([0.85, 0.05, 0.05, 0.35])

    extent_cross = [time.min(), time.max(), 0, max(cross_period)]
    extent_corr = [time.min(), time.max(), 0, max(cor_period)]
    im1 = NonUniformImage(ax1, interpolation='bilinear', extent=extent_cross)
    im1.set_data(time, cross_period, cross_power)
    ax1.images.append(im1)
    ax1.contour(time, cross_period, cross_sig, [-99, 1], colors='k', linewidths=2, extent=extent_cross)
    ax1.fill(np.concatenate([time, time[-1:] + dt, time[-1:] + dt, time[:1] - dt, time[:1] - dt]), np.concatenate([cross_coi, [1e-9], cross_period[-1:], cross_period[-1:], [1e-9]]), 'k', alpha=0.3, hatch='x')
    ax1.set_title('Cross-Wavelet')
    ax1.quiver(time[::3], cross_period[::3], u[::3, ::3], v[::3, ::3], units='width', angles='uv', pivot='mid', linewidth=1, edgecolor='k', headwidth=10, headlength=10, headaxislength=5, minshaft=2, minlength=5)
    fig.colorbar(im1, cax=cbar_ax)

    im2 = NonUniformImage(ax2, interpolation='bilinear', extent=extent_corr)
    im2.set_data(time, cor_period, WCT)
    ax2.images.append(im2)
    ax2.contour(time, cor_period, cor_sig, [-99, 1], colors='k', linewidths=2, extent=extent_corr)
    ax2.fill(np.concatenate([time, time[-1:] + dt, time[-1:] + dt, time[:1] - dt, time[:1] - dt]), np.concatenate([corr_coi, [1e-9], cor_period[-1:], cor_period[-1:], [1e-9]]),'k', alpha=0.3, hatch='x')
    ax2.set_title('Cross-Correlation')
    ax2.quiver(time[::3], cor_period[::3], u[::3, ::3], v[::3, ::3], units='height', angles='uv', pivot='mid', linewidth=1, edgecolor='k', headwidth=10, headlength=10, headaxislength=5, minshaft=2, minlength=5)
    #ax2.set_ylim(2, 35)
    #ax2.set_xlim(max(t1.min(), t2.min()), min(t1.max(), t2.max()))
    fig.colorbar(im2, cax=cbar_ax_1)

    plt.draw()
    plt.show()

    sys.exit(0)


if __name__ == "__main__":
    main()
