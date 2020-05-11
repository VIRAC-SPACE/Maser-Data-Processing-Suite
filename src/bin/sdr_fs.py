#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
SDR output data processing tool
"""

import sys
import scipy
import numpy as np
from PyQt5.QtWidgets import QWidget, QApplication
from PyQt5.QtGui import QIcon


def dopler(observed_frequency, velocity_receiver, base_frequency):
    """

    :param observed_frequency: observed frequency
    :param velocity_receiver: velocity of receiver
    :param base_frequency: laboratory determined frequency
    :return: velocity of source
    """
    speed_of_light = scipy.constants.speed_of_light
    velocity_source = (-((observed_frequency / base_frequency) - 1) *
                       speed_of_light + (velocity_receiver * 1000)) / 1000
    return velocity_source


def signal_to_noise_ratio(frequency, amplitude, cuts):
    """

    :param frequency: frequency
    :param amplitude: amplitude
    :param cuts: signal region
    :return: signal to noise ratio
    """
    cuts_index = list()
    cuts_index.append(0)

    for cut in cuts:
        cuts_index.append((np.abs(frequency - float(cut[0]))).argmin())
        cuts_index.append((np.abs(frequency - float(cut[1]))).argmin())

    cuts_index.append(-1)

    y_array = list()

    i = 0
    j = 1

    while i != len(cuts_index):
        y_array.append(amplitude[cuts_index[i]: cuts_index[j]])
        i = i + 2
        j = j + 2

    non_signal_amplitude = list()

    for point in y_array:
        for point_one in point:
            non_signal_amplitude.append(point_one)

    non_signal_amplitude = np.array(non_signal_amplitude)

    std = np.std(non_signal_amplitude)
    max_value = np.max(amplitude)

    ston = max_value / (std * 3)
    return ston


def frequency_shifting(p_sig_left, p_sig_right, p_ref_left, p_ref_right, p_sig_on_left,
                       p_sig_on_right, p_ref_on_left, p_ref_on_right, frequency_a, logs, pair):
    """

    :param p_sig_left: p_sig_left
    :param p_sig_right: p_sig_right
    :param p_ref_left: p_ref_left
    :param p_ref_right: p_ref_right
    :param p_sig_on_left: p_sig_on_left
    :param p_sig_on_right: p_sig_on_right
    :param p_ref_on_left: p_ref_on_left
    :param p_ref_on_right: p_ref_on_right
    :param frequency_a: frequency_a
    :param logs: logs
    :param pair: pair
    :return: frequency_shifting
    """
    df_div = float(logs["header"]["df_div,df"][0])
    band_with = float(logs["header"]["Fs,Ns,RBW"][0])
    f_shift = band_with / df_div

    l_spec = len(frequency_a)
    f_step = (frequency_a[l_spec - 1] - frequency_a[0]) / (l_spec - 1)
    n_shift = int(np.rint(f_shift / f_step))
    avg_interval = 0.5  # inner 50%
    s_i = int(l_spec / 2 - l_spec * avg_interval / 2)
    e_i = int(l_spec / 2 + l_spec * avg_interval / 2)

    tsys_off_1_left = float(logs["header"]["Tcal"][0]) * (
                (p_ref_on_left + p_ref_left) -
                np.mean(p_ref_on_left[s_i:e_i] - p_ref_left[s_i:e_i])) \
                / (2 * np.mean(p_ref_on_left[s_i:e_i] - p_ref_left[s_i:e_i]))
    tsys_off_2_left = float(logs["header"]["Tcal"][1]) * (
                (p_sig_on_left + p_sig_left) -
                np.mean(p_sig_on_left[s_i:e_i] -p_sig_left[s_i:e_i])) \
                / (2 * np.mean(p_sig_on_left[s_i:e_i] - p_sig_left[s_i:e_i]))

    tsys_off_1_right = float(logs["header"]["Tcal"][0]) * (
                (p_ref_on_right + p_ref_right) -
                np.mean(p_ref_on_right[s_i:e_i] - p_ref_right[s_i:e_i])) \
                / (2 * np.mean(p_ref_on_right[s_i:e_i] - p_ref_right[s_i:e_i]))
    tsys_off_2_right = float(logs["header"]["Tcal"][1]) * (
                (p_sig_on_right + p_sig_right) -
                np.mean(p_sig_on_right[s_i:e_i] - p_sig_right[s_i:e_i])) \
                / (2 * np.mean(p_sig_on_right[s_i:e_i] - p_sig_right[s_i:e_i]))

    ta_1_caloff_left = tsys_off_1_left * (p_sig_left - p_ref_left) / p_ref_left  # non-cal phase
    ta_1_caloff_right = tsys_off_1_right * \
                        (p_sig_right - p_ref_right) / p_ref_right  # non-cal phase

    ta_1_calon_left = (tsys_off_1_left + float(logs["header"]["Tcal"][0])) * \
                      (p_sig_on_left - p_ref_on_left) / p_ref_on_left  # cal phase
    ta_1_calon_right = (tsys_off_1_right + float(logs["header"]["Tcal"][1])) * \
                       (p_sig_on_right - p_ref_on_right) / p_ref_on_right  # cal phase

    ta_sig_left = (ta_1_caloff_left + ta_1_calon_left) / 2
    ta_sig_right = (ta_1_caloff_right + ta_1_calon_right) / 2

    ta_2_caloff_left = tsys_off_2_left * (p_ref_left - p_sig_left) / p_sig_left  # non-cal phase
    ta_2_caloff_right = tsys_off_2_right * \
                        (p_ref_right - p_sig_right) / p_sig_right  # non-cal phase

    ta_2_calon_left = (tsys_off_2_left + float(logs["header"]["Tcal"][0])) *\
                      (p_ref_on_left - p_sig_on_left) / p_sig_on_left  # cal phase
    ta_2_calon_right = (tsys_off_2_right + float(logs["header"]["Tcal"][1])) * \
                      (p_ref_on_right - p_sig_on_right) / p_sig_on_right  # cal phase

    ta_ref_left = (ta_2_caloff_left + ta_2_calon_left) / 2
    ta_ref_right = (ta_2_caloff_right + ta_2_calon_right) / 2

    ta_sig_left = np.roll(ta_sig_left, +n_shift)
    ta_sig_right = np.roll(ta_sig_right, +n_shift)

    ta_ref_left = np.roll(ta_ref_left, -n_shift)
    ta_ref_right = np.roll(ta_ref_right, -n_shift)

    ta_left = (ta_sig_left + ta_ref_left) / 2
    ta_right = (ta_sig_right + ta_ref_right) / 2

    tsys_r_left = np.mean(tsys_off_1_left[s_i:e_i])
    tsys_r_right = np.mean(tsys_off_1_right[s_i:e_i])

    tsys_s_left = np.mean(tsys_off_2_left[s_i:e_i])
    tsys_s_right = np.mean(tsys_off_2_right[s_i:e_i])

    elvation = (float(logs[pair[0][0]]["AzEl"][1]) + float(logs[pair[0][1]]["AzEl"][1]) + float(
        logs[pair[1][0]]["AzEl"][1]) + float(logs[pair[1][1]]["AzEl"][1])) / 4

    g_el = logs["header"]["Elev_poly"]
    g_el = [float(gel) for gel in g_el]
    g_e_ltmp = [0, 0, 0]
    g_e_ltmp[0] = g_el[2]
    g_e_ltmp[1] = g_el[1]
    g_e_ltmp[2] = g_el[0]
    g_el = g_e_ltmp

    sf_left = ta_left / ((float(logs["header"]["DPFU"][0])) * np.polyval(g_el, elvation))
    sf_right = ta_right / ((float(logs["header"]["DPFU"][1])) * np.polyval(g_el, elvation))

    return sf_left[s_i:e_i], sf_right[s_i:e_i], frequency_a[s_i:e_i], \
           tsys_r_left, tsys_r_right, tsys_s_left, tsys_s_right


class Analyzer(QWidget):
    """
    GUI application
    """
    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon('../../viraclogo.png'))
        self.setWindowTitle("SDR")


def main():
    """

    :return: None
    """
    q_app = QApplication(sys.argv)
    aplication = Analyzer()
    aplication.show()
    sys.exit(q_app.exec_())


if __name__ == "__main__":
    main()
