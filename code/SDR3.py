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

    Ns = 2048  # 32768/16;
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


def frequency_shifting(p_sig_left, p_sig_right, p_ref_left, p_ref_right, p_sig_on_left, p_sig_on_right, p_ref_on_left, p_ref_on_right, frequencyA, scan, logs):
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

    pair = ((str(scan) + "r0", str(scan) + "s0"), (str(scan) + "r1", str(scan) + "s1"))
    El = (float(logs[pair[0][0]]["AzEl"][1]) + float(logs[pair[0][1]]["AzEl"][1]) + float(logs[pair[1][0]]["AzEl"][1]) + float(logs[pair[1][1]]["AzEl"][1])) / 4
    G_El = logs["header"]["Elev_poly"]
    G_El = [float(gel) for gel in G_El]
    G_ELtmp = [0, 0, 0]
    G_ELtmp[0] = G_El[2]
    G_ELtmp[1] = G_El[1]
    G_ELtmp[2] = G_El[0]
    G_El = G_ELtmp

    Sf_left = Ta_left / (((-1) * float(logs["header"]["DPFU"][0])) * np.polyval(G_El, El))
    Sf_right = Ta_right / (((-1) * float(logs["header"]["DPFU"][1])) * np.polyval(G_El, El))

    si = si
    ei = ei

    return (Sf_left[si:ei], Sf_right[si:ei], frequencyA[si:ei])


def main():
    dir = "/home/janis/Documents/cepa_f6668_95/"
    scans = [1,2,3,4,5]
    logs = LogReaderFactory.getLogReader(LogTypes.SDR, "/home/janis/Downloads/cepa_f6668_95.log", "/home/janis/").getLogs()

    Sf_lefts = []
    Sf_rights = []

    Sf_lefts2 = []
    Sf_rights2 = []

    for scan in scans:
        # process data file

        file1 = dir + "cepa_f6668_95_no00" + str(scan) + "r0.dat"
        file2 = dir + "cepa_f6668_95_no00" + str(scan) + "s0.dat"
        file3 = dir + "cepa_f6668_95_no00" + str(scan) + "r1.dat"
        file4 = dir + "cepa_f6668_95_no00" + str(scan) + "s1.dat"

        frequencyA = read_dat(file1)[0]  # r0
        p_sig_left = read_dat(file1)[1]  # r0
        p_sig_right = read_dat(file1)[2]  # r0

        frequencyB = read_dat(file2)[0]  # s0
        p_ref_left = read_dat(file2)[1]  # s0
        p_ref_right = read_dat(file2)[2]  # s0

        frequencyC = read_dat(file3)[0]  # r1
        p_sig_on_left = read_dat(file3)[1]  # r1
        p_sig_on_right = read_dat(file3)[2]  # r1

        frequencyD = read_dat(file4)[0]  # s1
        p_ref_on_left = read_dat(file4)[1]  # s1
        p_ref_on_right = read_dat(file4)[2]  # s1

        p_sig_left = np.fft.fftshift(p_sig_left)  # r0
        p_sig_right = np.fft.fftshift(p_sig_right)  # r0
        p_ref_left = np.fft.fftshift(p_ref_left)  # s0
        p_ref_right = np.fft.fftshift(p_ref_right)  # s0
        p_sig_on_left = np.fft.fftshift(p_sig_on_left)  # r1
        p_sig_on_right = np.fft.fftshift(p_sig_on_right)  # r1
        p_ref_on_left = np.fft.fftshift(p_ref_on_left)  # s1
        p_ref_on_right = np.fft.fftshift(p_ref_on_right)  # s1

        Sf_left, Sf_right, frequencyA1 = frequency_shifting(p_sig_left, p_sig_right, p_ref_left, p_ref_right, p_sig_on_left, p_sig_on_right, p_ref_on_left, p_ref_on_right, frequencyA, scan, logs)
        Sf_lefts.append(Sf_left)
        Sf_rights.append(Sf_right)

        # process raw files
        file1 = dir + "cepa_f6668_95_no00" + str(scan) + "r0_ch1.raw"
        file2 = dir + "cepa_f6668_95_no00" + str(scan) + "s0_ch1.raw"
        file3 = dir + "cepa_f6668_95_no00" + str(scan) + "r1_ch1.raw"
        file4 = dir + "cepa_f6668_95_no00" + str(scan) + "s1_ch1.raw"

        file5 = dir + "cepa_f6668_95_no00" + str(scan) + "r0_ch2.raw"
        file6 = dir + "cepa_f6668_95_no00" + str(scan) + "s0_ch2.raw"
        file7 = dir + "cepa_f6668_95_no00" + str(scan) + "r1_ch2.raw"
        file8 = dir + "cepa_f6668_95_no00" + str(scan) + "s1_ch2.raw"

        read_raw(file1)
        read_raw(file2)
        read_raw(file3)
        read_raw(file4)
        read_raw(file5)
        read_raw(file6)
        read_raw(file7)
        read_raw(file8)

        p_sig_left = read_raw(file1)
        p_ref_left = read_raw(file2)
        p_sig_on_left = read_raw(file3)
        p_ref_on_left = read_raw(file4)

        p_sig_right = read_raw(file5)
        p_ref_right = read_raw(file6)
        p_sig_on_right = read_raw(file7)
        p_ref_on_right = read_raw(file8)

        print(len(p_sig_left), len(frequencyA))

        freq = []

        start = frequencyA[0]
        for f in range(0, 4096):
            freq.append(start)
            start += 0.000381

        freq = np.array(freq)

        Sf_left2, Sf_right2, frequencyA2 = frequency_shifting(p_sig_left, p_sig_right, p_ref_left, p_ref_right,p_sig_on_left, p_sig_on_right, p_ref_on_left, p_ref_on_right, freq, scan, logs)
        Sf_lefts2.append(Sf_left2)
        Sf_rights2.append(Sf_right2)

    # data files
    Sf_lefts_np = np.zeros(len(Sf_lefts[0]))
    for s in  Sf_lefts:
        Sf_lefts_np = Sf_lefts_np + np.array(s)

    Sf_rights_np = np.zeros(len(Sf_rights[0]))
    for s in  Sf_rights:
        Sf_rights_np = Sf_rights_np + np.array(s)

    Sf_lefts_np = Sf_lefts_np/len(Sf_lefts_np)
    Sf_rights_np = Sf_rights_np/len(Sf_rights_np)

    # raw files
    Sf_lefts_np2 = np.zeros(len(Sf_lefts2[0]))
    for s in Sf_lefts2:
        Sf_lefts_np2 = Sf_lefts_np2 + np.array(s)


    Sf_rights_np2 = np.zeros(len(Sf_rights2[0]))

    for s in Sf_rights2:
        Sf_rights_np2 = Sf_rights_np2 + np.array(s)

    Sf_lefts_np2 = Sf_lefts_np2 / len(Sf_lefts_np2)
    Sf_rights_np2 = Sf_rights_np2 / len(Sf_rights_np2)

    plt.subplot(2,2,1)
    plt.plot(frequencyA1, Sf_lefts_np)

    plt.subplot(2, 2, 2)
    plt.plot(frequencyA1, Sf_lefts_np)

    plt.subplot(2, 2, 3)
    plt.plot(Sf_lefts_np2)

    plt.subplot(2, 2, 4)
    plt.plot(Sf_lefts_np2)

    plt.show()

    sys.exit(0)


if __name__ == "__main__":
    main()