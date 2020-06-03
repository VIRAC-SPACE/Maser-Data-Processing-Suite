#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
SDR output data processing tool
"""

import sys
import os
import re
import argparse
import datetime
from functools import reduce
import scipy.constants
import numpy as np
from astropy.time import Time
import h5py
from PyQt5.QtWidgets import QWidget, QApplication, QDesktopWidget, QGridLayout, QPushButton
from PyQt5.QtGui import QIcon
from ExperimentsLogReader.experimentsLogReader import LogReaderFactory, LogTypes
from parsers.configparser_ import ConfigParser
from utils.vlsr import lsr
from utils.help import find_nearest_index
from utils.ploting_qt5 import Plot


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Creates input file for plotting tool. ''',
                                     epilog="""PRE PLOTTER.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("line", help="frequency", type=str)
    parser.add_argument("iteration_number", help="iteration number ", type=int)
    parser.add_argument("log_file", help="Experiment log file name", type=str)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str,
                        default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    return args


def get_args(key):
    """

    :param key: argument key
    :return: to script passed argument value
    """
    return str(parse_arguments().__dict__[key])


def get_configs(section, key):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :return: configuration file section key value
    """
    config_file_path = get_args("config")
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


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

    tsys_off_1_left = float(logs["header"]["Tcal"][0]) * ((p_ref_on_left + p_ref_left) -
                                                          np.mean(p_ref_on_left[s_i:e_i] -
                                                                  p_ref_left[s_i:e_i])) / \
                      (2 * np.mean(p_ref_on_left[s_i:e_i] - p_ref_left[s_i:e_i]))
    tsys_off_2_left = float(logs["header"]["Tcal"][1]) * ((p_sig_on_left + p_sig_left)
                                                          - np.mean(p_sig_on_left[s_i:e_i]
                                                                    - p_sig_left[s_i:e_i])) / \
                      (2 * np.mean(p_sig_on_left[s_i:e_i] - p_sig_left[s_i:e_i]))

    tsys_off_1_right = float(logs["header"]["Tcal"][0]) * ((p_ref_on_right + p_ref_right) -
                                                           np.mean(p_ref_on_right[s_i:e_i]
                                                                   - p_ref_right[s_i:e_i])) / \
                       (2 * np.mean(p_ref_on_right[s_i:e_i] - p_ref_right[s_i:e_i]))

    tsys_off_2_right = float(logs["header"]["Tcal"][1]) * ((p_sig_on_right + p_sig_right) -
                                                           np.mean(p_sig_on_right[s_i:e_i]
                                                                   - p_sig_right[s_i:e_i])) / \
                       (2 * np.mean(p_sig_on_right[s_i:e_i] - p_sig_right[s_i:e_i]))

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

    ta_2_calon_left = (tsys_off_2_left + float(logs["header"]["Tcal"][0])) * \
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


def get_scan_name(data_file_name):
    """

    :param data_file_name: data file name
    :return: scan name for data file
    """
    return data_file_name.split(".")[0].split("_")[4][2:len(data_file_name.split(".")
                                                            [0].split("_")[4])].lstrip('0')


def get_data(data_file_name):
    """

    :param data_file_name: data file name
    :return: frequency, polarization left, polarization right
    """
    frequency = np.loadtxt(data_file_name, usecols=(0,), unpack=True)
    polarization_left = np.loadtxt(data_file_name, usecols=(1,), unpack=True)
    polarization_right = np.loadtxt(data_file_name, usecols=(2,), unpack=True)
    return frequency, polarization_left, polarization_right


class Analyzer(QWidget):
    """
    GUI application
    """

    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon('viraclogo.png'))
        self.setWindowTitle("SDR")
        self.center()
        self.index = 0
        self.cuts = get_configs('cuts', get_args("source") + "_" + get_args("line")).split(";")
        self.cuts = [c.split(",") for c in self.cuts]
        self.sf_left = list()
        self.sf_right = list()
        self.tsys_r_left_list = list()
        self.tsys_r_right_list = list()
        self.tsys_s_left_list = list()
        self.tsys_s_right_list = list()
        self.ston_list_left = list()
        self.ston_list_right = list()
        self.ston_list_avg = list()
        self.logs = LogReaderFactory.getLogReader(LogTypes.SDR,
                                                  get_configs("paths", "logPath") + "SDR/" +
                                                  get_args("log_file"),
                                                  get_configs("paths", "prettyLogsPath") +
                                                  get_args("source") + "_" +
                                                  get_args("iteration_number")).getLogs()
        self.station = self.logs["header"]["station,id"]
        self.data_dir = get_configs("paths", "dataFilePath") + \
                        get_args("source") + "_f" + get_args("line") + "_" + \
                        self.station[1] + "_" + get_args("iteration_number") + "/"
        self.data_files = os.listdir(self.data_dir)
        self.scan_pairs = self.create_scan_pairs()
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.setSpacing(10)

        self.plot_velocity__left = None
        self.plot_velocity__right = None
        self.plot_tsys = None
        self.plot_ston = None
        self.plot_start__left_a = None
        self.plot_start__right_b = None
        self.x = None
        self.total__left = None
        self.total__right = None

        self.__UI__()

    def center(self):
        """

        :return: None
        """
        frame_geometry = self.frameGeometry()
        centre_position = QDesktopWidget().availableGeometry().center()
        frame_geometry.moveCenter(centre_position)
        self.move(frame_geometry.topLeft())

    def __UI__(self):
        """

        :return: None
        """
        if self.index != len(self.scan_pairs) - 1:  # cheking if there is not one pair
            self.next_pair_button = QPushButton("Next pair", self)
            self.next_pair_button.clicked.connect(self.next_pair)
            self.grid.addWidget(self.next_pair_button, 4, 2)

        self.skip_all_button = QPushButton("Skip to end", self)
        self.skip_all_button.clicked.connect(self.skip_all)
        self.grid.addWidget(self.skip_all_button, 5, 2)

        self.plot_pair(self.index)

    def create_scan_pairs(self):
        """

        :return: None
        """
        scan_names = [get_scan_name(file) for file in self.data_files]
        scans_numbers = list(set([int(re.findall("[0-9]+", s)[0]) for s in scan_names]))
        scans_numbers = sorted(scans_numbers)
        scan_pairs = []

        for scan in scans_numbers:
            scan_pairs.append(
                ((str(scan) + "r" + "0", str(scan) + "s" + "0"),
                 (str(scan) + "r" + "1", str(scan) + "s" + "1")))

        return scan_pairs

    def next_pair(self):
        """

        :return: None
        """
        if self.index == len(self.scan_pairs) - 1:
            pass

        else:
            self.plot_start__left_a.removePolt()
            self.plot_start__right_b.removePolt()
            self.index = self.index + 1
            self.plot_pair(self.index)

    def skip_all(self):
        """

        :return: None
        """
        self.index += 1
        while self.index < len(self.scan_pairs):
            pair = self.scan_pairs[self.index]

            file1 = self.data_dir + self.get_data_file_for_scan(pair[0][1])  # s0
            file2 = self.data_dir + self.get_data_file_for_scan(pair[0][0])  # r0
            file3 = self.data_dir + self.get_data_file_for_scan(pair[1][1])  # s1
            file4 = self.data_dir + self.get_data_file_for_scan(pair[1][0])  # r1

            print("data files", file1, file2, file3, file4)

            frequency_a = get_data(file1)[0]  # r0
            p_sig_left = get_data(file1)[1]  # r0
            p_sig_right = get_data(file1)[2]  # r0

            p_ref_left = get_data(file2)[1]  # s0
            p_ref_right = get_data(file2)[2]  # s0

            p_sig_on_left = get_data(file3)[1]  # r1
            p_sig_on_right = get_data(file3)[2]  # r1

            p_ref_on_left = get_data(file4)[1]  # s1
            p_ref_on_right = get_data(file4)[2]  # s1

            # fft shift
            p_sig_left = np.fft.fftshift(p_sig_left)  # r0
            p_sig_right = np.fft.fftshift(p_sig_right)  # r0
            p_ref_left = np.fft.fftshift(p_ref_left)  # s0
            p_ref_right = np.fft.fftshift(p_ref_right)  # s0
            p_sig_on_left = np.fft.fftshift(p_sig_on_left)  # r1
            p_sig_on_right = np.fft.fftshift(p_sig_on_right)  # r1
            p_ref_on_left = np.fft.fftshift(p_ref_on_left)  # s1
            p_ref_on_right = np.fft.fftshift(p_ref_on_right)  # s1

            sf_left, sf_right, frequency_a1, tsys_r_left, \
            tsys_r_right, tsys_s_left, tsys_s_right = frequency_shifting(
                p_sig_left, p_sig_right, p_ref_left, p_ref_right,
                p_sig_on_left, p_sig_on_right, p_ref_on_left,
                p_ref_on_right, frequency_a, self.logs, pair)

            self.sf_left.append(sf_left)
            self.sf_right.append(sf_right)
            self.x = frequency_a1

            self.tsys_r_left_list.append(tsys_r_left)
            self.tsys_r_right_list.append(tsys_r_right)
            self.tsys_s_left_list.append(tsys_s_left)
            self.tsys_s_right_list.append(tsys_s_right)

            ston_left = signal_to_noise_ratio(self.x, self.sf_left, self.cuts)
            ston_right = signal_to_noise_ratio(self.x, self.sf_right, self.cuts)
            stone_avg = signal_to_noise_ratio(self.x, ((np.array(self.sf_left) +
                                                        np.array(self.sf_right)) / 2), self.cuts)

            self.ston_list_left.append(ston_left)
            self.ston_list_right.append(ston_right)
            self.ston_list_avg.append(stone_avg)

            if self.index == len(self.scan_pairs) - 1:
                self.next_pair_button.setText('Move to total results')
                self.next_pair_button.clicked.connect(self.plot_total_results)

            self.index += 1

        self.plot_total_results()

    def plot_total_results(self):
        """

        :return: None
        """
        self.grid.removeWidget(self.plot_start__left_a)
        self.grid.removeWidget(self.plot_start__right_b)
        self.grid.removeWidget(self.total__left)
        self.grid.removeWidget(self.total__right)

        self.plot_start__left_a.hide()
        self.plot_start__right_b.hide()
        self.total__left.hide()
        self.total__right.hide()

        self.plot_start__left_a.close()
        self.plot_start__right_b.close()
        self.total__left.close()
        self.total__right.close()

        self.plot_start__left_a.removePolt()
        self.plot_start__right_b.removePolt()
        self.total__left.removePolt()
        self.total__right.removePolt()

        del self.plot_start__left_a
        del self.plot_start__right_b
        del self.total__left
        del self.total__right

        self.grid.removeWidget(self.next_pair_button)
        self.next_pair_button.hide()
        self.next_pair_button.close()
        del self.next_pair_button

        for i in reversed(range(self.grid.count())):
            self.grid.itemAt(i).widget().deleteLater()

        station = self.logs["header"]["station,id"][0]
        if station == "RT-32":
            station = "IRBENE"
        else:
            station = "IRBENE16"

        station_coordinates = get_configs("stations", station)
        station_coordinates = station_coordinates.replace(" ", "").split(",")
        x = np.float64(station_coordinates[0])
        y = np.float64(station_coordinates[1])
        z = np.float64(station_coordinates[2])

        velocity_max = []
        velocity_min = []
        velocity_list = []

        for p in range(0, len(self.scan_pairs)):
            scan_number = self.scan_pairs[p][0][0]
            scan_1 = self.logs[str(scan_number)]
            string_time = scan_1["date"].replace("T", " ")
            t = datetime.datetime.strptime(scan_1["date"], '%Y-%m-%dT%H:%M:%S')
            time = t.isoformat()
            date = Time(time, format='isot', scale='utc')

            print("station_coordinates", station_coordinates)

            source_cordinations = get_configs("sources", get_args("source")).split(",")
            source_cordinations = [sc.strip() for sc in source_cordinations]
            RA = source_cordinations[0]
            DEC = source_cordinations[1]

            ra = list()
            dec = list()
            ra.append(RA[0:2])
            ra.append(RA[2:4])
            ra.append(RA[4:len(RA)])

            if DEC[0] == "-":
                dec.append(DEC[0:3])
                dec.append(DEC[3:5])
                dec.append(DEC[5:len(DEC)])
            else:
                dec.append(DEC[0:2])
                dec.append(DEC[2:4])
                dec.append(DEC[4:len(DEC)])

            ra_str = ra[0] + "h" + ra[1] + "m" + ra[2] + "s"
            if int(dec[0]) > 0:
                dec_str = "+" + dec[0] + "d" + dec[1] + "m" + dec[2] + "s"
            else:
                dec_str = dec[0] + "d" + dec[1] + "m" + dec[2] + "s"

            print("Vel Total params", ra_str, dec_str, date, string_time, x, y, z)
            vel_total = lsr(ra_str, dec_str, date, string_time, x, y, z)
            print("vel_total", vel_total)

            line = get_configs('base_frequencies_SDR', "f" +
                               get_args("line")).replace(" ", "").split(",")
            line_f = float(line[0]) * (10 ** 9)
            line_s = line[1]
            specie = line_s

            print("specie", specie, "\n")
            local_oscillator = float(self.logs["header"]["f_obs,LO,IF"][1])
            velocities = dopler((self.x + local_oscillator) * (10 ** 6), vel_total, line_f)
            velocity_list.append(velocities)

            print("Velocity max value is", np.max(velocities),
                  "velocity min value is ", np.min(velocities))
            velocity_max.append(np.max(velocities))
            velocity_min.append(np.min(velocities))

        velocities_avg = []
        y__left_avg = []
        y__right_avg = []

        print("Velocity minimum from max value is", np.min(velocity_max),
              "velocity maximum for min value is ", np.max(velocity_min))
        left_cut = np.max(velocity_min)
        right_cut = np.min(velocity_max)

        for p in range(0, len(velocity_list)):
            index_left = find_nearest_index(velocity_list[p], left_cut)
            index_right = find_nearest_index(velocity_list[p], right_cut)
            print("index left ", index_left, "index right", index_right)

            y__left_avg.append(self.sf_left[p][index_right:index_left])
            y__right_avg.append(self.sf_right[p][index_right:index_left])
            velocities_avg.append(velocity_list[p][index_right:index_left])
            print(len(velocity_list[p][index_right:index_left]))

        max_points_count = np.max([len(m) for m in velocities_avg])
        print("max_points_count", max_points_count)
        for s in range(0, len(velocity_list)):

            if len(velocities_avg[s]) < max_points_count:
                velocities_avg[s] = np.append(velocities_avg[s], np.max(velocity_min))

            if len(y__left_avg[s]) < max_points_count:
                y__left_avg[s] = np.append(y__left_avg[s], 0)

            if len(y__right_avg[s]) < max_points_count:
                y__right_avg[s] = np.append(y__right_avg[s], 0)

        velocities_avg = reduce(lambda x, y: x + y, velocities_avg)
        y__left_avg = reduce(lambda x, y: x + y, y__left_avg)
        y__right_avg = reduce(lambda x, y: x + y, y__right_avg)

        velocities_avg = velocities_avg / len(self.scan_pairs)
        y__left_avg = y__left_avg / len(self.scan_pairs)
        y__right_avg = y__right_avg / len(self.scan_pairs)

        self.plot_velocity__left = Plot()
        self.plot_velocity__left.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)',
                                           'Flux density (Jy)',
                                           "Left Polarization", (1, 0), "linear")
        self.plot_velocity__left.plot(velocities_avg, y__left_avg, 'b')

        self.plot_velocity__right = Plot()
        self.plot_velocity__right.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)',
                                            'Flux density (Jy)',
                                            "Right Polarization", (1, 1), "linear")
        self.plot_velocity__right.plot(velocities_avg, y__right_avg, 'b')

        self.plot_tsys = Plot()
        self.plot_tsys.creatPlot(self.grid, 'Time', 'System temperature',
                                 "System temperature in time", (3, 0), "linear")

        time = list(set([int(t.split("_")[-1].split(".")[0]
                             [2:len(t.split("_")[-1].split(".")[0]) - 2])
                         for t in self.data_files]))
        self.plot_tsys.plot(time, self.tsys_r_left_list, '*b', label="Tsys_r_left")
        self.plot_tsys.plot(time, self.tsys_r_right_list, '*r', label="Tsys_r_right")
        self.plot_tsys.plot(time, self.tsys_s_left_list, '*g', label="Tsys_s_left")
        self.plot_tsys.plot(time, self.tsys_s_right_list, '*y', label="Tsys_s_right")

        # source_day_Month_year_houre:minute:righte_station_iteration.dat
        day = scan_1["date"].split("-")[2][0:2]
        month = scan_1["date"].split("-")[1]
        months = {"Jan": "1", "Feb": "2", "Mar": "3", "Apr": "4", "May": "5",
                  "Jun": "6", "Jul": "7", "Aug": "8", "Sep": "9", "Oct": "10",
                  "Nov": "11", "Dec": "12"}
        month = list(months.keys())[int(month) - 1]
        year = scan_1["date"].split("-")[0]
        hour = scan_1["date"].split("T")[1].split(":")[0]
        minute = scan_1["date"].split("T")[1].split(":")[1]
        right = scan_1["date"].split("T")[1].split(":")[2]

        if not os.path.exists(get_configs("paths", "outputFilePath") + "/" + get_args("line")):
            os.makedirs(get_configs("paths", "outputFilePath") + "/" + get_args("line"))

        result_file_name = get_configs("paths", "outputFilePath") + "/" + \
                           get_args("line") + "/" + \
                           get_args("source") + "_" + day + "_" + month + "_" + year + "_" + \
                           hour + ":" + minute + ":" + right + "_" + station + "_" + \
                           get_args("iteration_number") + ".h5"

        result_file = h5py.File(result_file_name, "w")

        sys_temp_out = np.transpose(np.array([time, self.tsys_r_left_list,
                                              self.tsys_r_right_list,
                                              self.tsys_s_left_list, self.tsys_s_right_list]))

        result_file.create_dataset("system_temperature", data=sys_temp_out)

        self.ston_list_left = [value for value in self.ston_list_left if str(value) != 'nan']
        self.ston_list_right = [value for value in self.ston_list_right if str(value) != 'nan']
        self.ston_list_avg = [value for value in self.ston_list_avg if str(value) != 'nan']

        time = list(time)
        while len(time) != len(self.ston_list_left):
            time.pop()

        self.plot_ston = Plot()
        self.plot_ston.creatPlot(self.grid, 'Pair', 'Ratio', "Signal to Noise", (3, 1), "linear")
        self.plot_ston.plot(time, self.ston_list_left, '*r', label="left Polarization")
        self.plot_ston.plot(time, self.ston_list_right, 'og', label="fight Polarization")
        self.plot_ston.plot(time, self.ston_list_avg, 'vb', label="AVG Polarization")

        print("Average signal to noise for left polarization", np.mean(self.ston_list_left))
        print("Average signal to noise for right polarization", np.mean(self.ston_list_right))
        print("Average signal to noise for average polarization", np.mean(self.ston_list_avg))

        self.grid.addWidget(self.plot_velocity__left, 0, 0)
        self.grid.addWidget(self.plot_velocity__right, 0, 1)
        self.grid.addWidget(self.plot_tsys, 2, 0)
        self.grid.addWidget(self.plot_ston, 2, 1)

        total_results = np.transpose(np.array([np.transpose(velocities_avg),
                                               np.transpose(y__left_avg),
                                               np.transpose(y__right_avg)]))

        result_file.create_dataset("amplitude", data=total_results)
        print("specie", specie)
        specie = [specie.encode("ascii", "ignore")]
        result_file.create_dataset("specie", (len(specie),1),'S10', specie)
        result_file.close()

    def get_data_file_for_scan(self, scan_name):
        """

        :param scan_name: scan name
        :return: file name for scan
        """
        file_name = ""

        for file in self.data_files:
            if get_scan_name(file) == scan_name:
                file_name = file
                break

        return file_name

    def plot_pair(self, index):
        """

        :param index: index of scans
        :return: None
        """
        pair = self.scan_pairs[index]
        file1 = self.data_dir + self.get_data_file_for_scan(pair[0][1])  # s0
        file2 = self.data_dir + self.get_data_file_for_scan(pair[0][0])  # r0
        file3 = self.data_dir + self.get_data_file_for_scan(pair[1][1])  # s1
        file4 = self.data_dir + self.get_data_file_for_scan(pair[1][0])  # r1

        print("data files", file1, file2, file3, file4)

        frequency_a, p_sig_left, p_sig_right = get_data(file1)  # s0
        frequency_b, p_ref_left, p_ref_right = get_data(file2)  # r0
        frequency_c, p_sig_on_left, p_sig_on_right = get_data(file3)  # s1
        frequency_d, p_ref_on_left, p_ref_on_right = get_data(file4)  # r1

        # fft shift
        p_sig_left = np.fft.fftshift(p_sig_left)  # s0
        p_sig_right = np.fft.fftshift(p_sig_right)  # s0
        p_ref_left = np.fft.fftshift(p_ref_left)  # r0
        p_ref_right = np.fft.fftshift(p_ref_right)  # r0
        p_sig_on_left = np.fft.fftshift(p_sig_on_left)  # s1
        p_sig_on_right = np.fft.fftshift(p_sig_on_right)  # s1
        p_ref_on_left = np.fft.fftshift(p_ref_on_left)  # r1
        p_ref_on_right = np.fft.fftshift(p_ref_on_right)  # r1

        sf_left, sf_right, frequency_a1, tsys_r_left, tsys_r_right, tsys_s_left, tsys_s_right = \
            frequency_shifting(p_sig_left, p_sig_right, p_ref_left, p_ref_right, p_sig_on_left,
                               p_sig_on_right, p_ref_on_left,
                               p_ref_on_right, frequency_a, self.logs, pair)

        self.sf_left.append(sf_left)
        self.sf_right.append(sf_right)
        self.x = frequency_a1

        self.tsys_r_left_list.append(tsys_r_left)
        self.tsys_r_right_list.append(tsys_r_right)
        self.tsys_s_left_list.append(tsys_s_left)
        self.tsys_s_right_list.append(tsys_s_right)

        ston_left = signal_to_noise_ratio(self.x, self.sf_left, self.cuts)
        ston_right = signal_to_noise_ratio(self.x, self.sf_right, self.cuts)
        stone_avg = signal_to_noise_ratio(self.x,
                                          ((np.array(self.sf_left) + np.array(self.sf_right)) / 2),
                                          self.cuts)

        self.ston_list_left.append(ston_left)
        self.ston_list_right.append(ston_right)
        self.ston_list_avg.append(stone_avg)

        # plot1
        self.plot_start__left_a = Plot()
        self.plot_start__left_a.creatPlot(self.grid, 'Frequency Mhz',
                                          'Amplitude', "Left Polarization", (1, 0), "linear")
        self.plot_start__left_a.plot(frequency_a, p_sig_left, 'b', label=str(index + 1) + "s0")
        self.plot_start__left_a.plot(frequency_b, p_ref_left, 'g', label=str(index + 1) + "r0")
        self.plot_start__left_a.plot(frequency_c, p_sig_on_left, 'r', label=str(index + 1) + "s1")
        self.plot_start__left_a.plot(frequency_d, p_ref_on_left, 'y', label=str(index + 1) + "r1")
        self.grid.addWidget(self.plot_start__left_a, 0, 0)

        # plot2
        self.plot_start__right_b = Plot()
        self.plot_start__right_b.creatPlot(self.grid, 'Frequency Mhz',
                                           'Amplitude', "Right Polarization", (1, 1), "linear")
        self.plot_start__right_b.plot(frequency_a, p_sig_right, 'b', label=str(index + 1) + "s0")
        self.plot_start__right_b.plot(frequency_b, p_ref_right, 'g', label=str(index + 1) + "r0")
        self.plot_start__right_b.plot(frequency_c, p_sig_on_right, 'r', label=str(index + 1) + "s1")
        self.plot_start__right_b.plot(frequency_d, p_ref_on_right, 'y', label=str(index + 1) + "r1")
        self.grid.addWidget(self.plot_start__right_b, 0, 1)

        # plot3
        self.total__left = Plot()
        self.total__left.creatPlot(self.grid, 'Frequency Mhz',
                                   'Flux density (Jy)', "", (4, 0), "linear")
        self.total__left.plot(self.x, sf_left, 'b', label=str(index + 1))
        self.grid.addWidget(self.total__left, 3, 0)

        # plot4
        self.total__right = Plot()
        self.total__right.creatPlot(self.grid,
                                    'Frequency Mhz', 'Flux density (Jy)', "", (4, 1), "linear")
        self.total__right.plot(self.x, sf_right, 'b', label=str(index + 1))
        self.grid.addWidget(self.total__right, 3, 1)

        if index == len(self.scan_pairs) - 1:
            self.next_pair_button.setText('Move to total results')
            self.next_pair_button.clicked.connect(self.plot_total_results)


def main():
    """

    :return: None
    """
    q_app = QApplication(sys.argv)
    application = Analyzer()
    application.show()
    sys.exit(q_app.exec_())


if __name__ == "__main__":
    main()
