#! /usr/bin/python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

from ExperimentsLogReader.experimentsLogReader import LogReaderFactory, LogTypes

from help import *


def read_raw(file_name):
    data = np.fromfile(file_name, np.dtype([('re', np.int16), ('im', np.int16)]))
    iq_data = data['re'] + data['im'] * 1j
    iq_data = iq_data * 1e-3

    Ns_tot = len(iq_data);
    #print("Ns_tot: %d samples" % (Ns_tot));

    # Tint = 50
    # Fs = 125e6
    # print("Expect: %d samples" % (int(Tint*Fs)));

    Ns = 4096  # 32768/16;
    Nint = int(Ns_tot / Ns);

    spec = np.zeros(Ns)
    ptr = 0;
    for i in range(Nint):
        spec_i = np.fft.fft(iq_data[ptr:(ptr + Ns)]);
        # spec = spec+10*nplog10(abs(spec_i))
        spec = spec + (spec_i.real * spec_i.real + spec_i.imag * spec_i.imag)
        ptr = ptr + Ns

    spec = np.fft.fftshift(spec)
    spec = spec / Ns;

    return spec


def read_dat(file_name):
    data = np.fromfile(file_name, dtype="float64", count=-1, sep=" ").reshape((file_len(file_name), 3))
    frequency = correctNumpyReadData(data[:, [0]])
    polarization_left = correctNumpyReadData(data[:, [1]])
    polarization_right = correctNumpyReadData(data[:, [2]])
    return (frequency, polarization_left, polarization_right)


def frequency_shifting(p_sig_left, p_sig_right, p_ref_left, p_ref_right, p_sig_on_left, p_sig_on_right, p_ref_on_left, p_ref_on_right, frequencyA, logs):
    df_div = float(logs["header"]["df_div,df"][0])
    BW = float(logs["header"]["Fs,Ns,RBW"][0])
    f_shift = BW / df_div

    l_spec = len(frequencyA)
    f_step = (frequencyA[l_spec - 1] - frequencyA[0]) / (l_spec - 1)
    n_shift = int(np.rint(f_shift / f_step))
    avg_interval = 0.5  # inner 50%
    si = int(l_spec / 2 - l_spec * avg_interval / 2)
    ei = int(l_spec / 2 + l_spec * avg_interval / 2)

    Tsys_off_1_left = float(logs["header"]["Tcal"][0]) * ((p_ref_on_left + p_ref_left) - np.mean(p_ref_on_left[si:ei] - p_ref_left[si:ei])) / (2 * np.mean(p_ref_on_left[si:ei] - p_ref_left[si:ei]))
    Tsys_off_2_left = float(logs["header"]["Tcal"][1]) * ((p_sig_on_left + p_sig_left) - np.mean(p_sig_on_left[si:ei] - p_sig_left[si:ei])) / (2 * np.mean(p_sig_on_left[si:ei] - p_sig_left[si:ei]))

    Tsys_off_1_right = float(logs["header"]["Tcal"][0]) * ((p_ref_on_right + p_ref_right) - np.mean(p_ref_on_right[si:ei] - p_ref_right[si:ei])) / (2 * np.mean(p_ref_on_right[si:ei] - p_ref_right[si:ei]))
    Tsys_off_2_right = float(logs["header"]["Tcal"][1]) * ((p_sig_on_right + p_sig_right) - np.mean(p_sig_on_right[si:ei] - p_sig_right[si:ei])) / (2 * np.mean(p_sig_on_right[si:ei] - p_sig_right[si:ei]))

    Ta_1_caloff_left = Tsys_off_1_left * (p_sig_left - p_ref_left) / p_ref_left  # non-cal phase
    Ta_1_caloff_right = Tsys_off_1_right * (p_sig_right - p_ref_right) / p_ref_right  # non-cal phase

    Ta_1_calon_left = (Tsys_off_1_left + float(logs["header"]["Tcal"][0])) * (p_sig_on_left - p_ref_on_left) / p_ref_on_left  # cal phase
    Ta_1_calon_right = (Tsys_off_1_right + float(logs["header"]["Tcal"][1])) * (p_sig_on_right - p_ref_on_right) / p_ref_on_right  # cal phase

    Ta_sig_left = (Ta_1_caloff_left + Ta_1_calon_left) / 2
    Ta_sig_right = (Ta_1_caloff_right + Ta_1_calon_right) / 2

    Ta_2_caloff_left = Tsys_off_2_left * (p_ref_left - p_sig_left) / p_sig_left  # non-cal phase
    Ta_2_caloff_right = Tsys_off_2_right * (p_ref_right - p_sig_right) / p_sig_right  # non-cal phase

    Ta_2_calon_left = (Tsys_off_2_left + float(logs["header"]["Tcal"][0])) * (p_ref_on_left - p_sig_on_left) / p_sig_on_left  # cal phase
    Ta_2_calon_right = (Tsys_off_2_right + float(logs["header"]["Tcal"][1])) * (p_ref_on_right - p_sig_on_right) / p_sig_on_right  # cal phase

    Ta_ref_left = (Ta_2_caloff_left + Ta_2_calon_left) / 2
    Ta_ref_right = (Ta_2_caloff_right + Ta_2_calon_right) / 2

    Ta_sig_left = np.roll(Ta_sig_left, +n_shift)
    Ta_sig_right = np.roll(Ta_sig_right, +n_shift)

    Ta_ref_left = np.roll(Ta_ref_left, -n_shift)
    Ta_ref_right = np.roll(Ta_ref_right, -n_shift)

    Ta_left = (Ta_sig_left + Ta_ref_left) / 2
    Ta_right = (Ta_sig_right + Ta_ref_right) / 2

    Tsys_r_left = np.mean(Tsys_off_1_left[si:ei])
    Tsys_r_right = np.mean(Tsys_off_1_right[si:ei])

    Tsys_s_left = np.mean(Tsys_off_2_left[si:ei])
    Tsys_s_right = np.mean(Tsys_off_2_right[si:ei])

    El = (float(logs[pair[0][0]]["AzEl"][1]) + float(logs[pair[0][1]]["AzEl"][1]) + float(logs[pair[1][0]]["AzEl"][1]) + float(logs[pair[1][1]]["AzEl"][1])) / 4

    G_El = logs["header"]["Elev_poly"]
    G_El = [float(gel) for gel in G_El]
    G_ELtmp = [0, 0, 0]
    G_ELtmp[0] = G_El[2]
    G_ELtmp[1] = G_El[1]
    G_ELtmp[2] = G_El[0]
    G_El = G_ELtmp

    Sf_left = Ta_left / ((float(logs["header"]["DPFU"][0])) * np.polyval(G_El, El))
    Sf_right = Ta_right / ((float(logs["header"]["DPFU"][1])) * np.polyval(G_El, El))

    return (Sf_left[si:ei], Sf_right[si:ei], frequencyA[si:ei])


def main():
    raw_file_dir = "/mnt/WORK/temp/"
    data_file_dir = 
    
    sys.exit(0)


if __name__ == "__main__":
    main()
