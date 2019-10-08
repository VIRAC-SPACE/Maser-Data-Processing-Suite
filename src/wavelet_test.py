#! /usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from obspy.imaging.cm import obspy_sequential
from obspy.signal.tf_misfit import cwt
import argparse

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
dt = (t0 - tmax)/len(time)
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

