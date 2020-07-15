#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
creates spectre movie
"""

import sys
import os
import argparse
import subprocess
import json
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QApplication, QWidget, QGridLayout, QDesktopWidget, QLabel, QToolButton, QPushButton
import h5py
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from utils.ploting_qt5 import Plot
from parsers.configparser_ import ConfigParser


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Monitoring spectr change in time. ''',
                                     epilog="""Monitor Spectr.""")
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("line", help="line", type=int)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="../config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    parser.add_argument("-i", "--index", help="Configuration cfg file", type=int, default=0)
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


class SpectreTime(QWidget):
    def __init__(self, source_name):
        super().__init__()
        self.setWindowIcon(QIcon('../viraclogo.png'))
        self.center()
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.source_name = source_name
        self.setWindowTitle(self.source_name)

        self.output_path = get_configs("paths", "outputFilePath") + "/" + get_args("line") + "/"
        output_files = [file for file in os.listdir(self.output_path) if file.startswith(self.source_name + "_")]
        self.index = int(get_args("index"))
        dates = sorted([file_name.split("_")[1] for file_name in output_files])
        self.sorted_file_names = []
        for f in range(0, len(dates)):
            for file_name in output_files:
                if file_name.split("_")[1] == dates[f]:
                    self.sorted_file_names.append(file_name)

        self.information_label = QLabel("")
        self.grid.addWidget(self.information_label, 0, 2)
        self.previous_button = QToolButton(arrowType=Qt.LeftArrow)
        self.next_button = QToolButton(arrowType=Qt.RightArrow)
        self.previous_button.clicked.connect(self.previous_spectre)
        self.next_button.clicked.connect(self.next_spectre)
        self.grid.addWidget(self.previous_button, 1, 0)
        self.grid.addWidget(self.next_button, 1, 2)
        create_movie_button = QPushButton('Create movie', self)
        create_movie_button.clicked.connect(self.create_movie)
        self.grid.addWidget(create_movie_button, 1, 3)
        self.x_lim = None
        self.y_lim = None
        self.previous_line = None
        self.first_plot = False
        result_file = get_configs("paths", "resultFilePath") + "/" + self.source_name + "_" + get_args("line") + ".json"
        with open(result_file) as result_data:
            self.results = json.load(result_data)

        self.spectre_plot = Plot()
        self.spectre_plot.creatPlot(self.grid, "Velocity (km sec$^{-1}$)", "Flux density (Jy)",
                                    get_configs("Full_source_name", self.source_name), (1, 1), "linear")
        self.spectre_plot.addZoomEvent(self.zoom_callback)
        self.plot()

    def plot(self):
        """

        :return: None
        """
        file_name = self.output_path + self.sorted_file_names[self.index]
        data = h5py.File(file_name, 'r')['amplitude_corrected'][()]
        xdata = data[:, 0]
        ydata = data[:, 3]

        modified_julian_days = file_name.split("_")[1]
        date = Time(modified_julian_days, format='mjd').strftime('%d %b %Y')
        iteration = file_name.split("_")[3]
        experiment_name = ".".join([file_name.split("/")[-1].split(".")[0],
                              file_name.split("/")[-1].split(".")[1]])

        type_of_observation = self.results[experiment_name]["type"]
        flag = self.results[experiment_name]["flag"]
        self.information_label.setText("Date " + date + "\n" + "Iteration " + str(iteration) +
                                       "\n" + "Modified Julian Days " + str(modified_julian_days) + "\n" +
                                       "Type " + type_of_observation + "\n" + "Flag " + str(flag))

        line = self.spectre_plot.plot(xdata, ydata, "-")

        if self.previous_line:
            self.previous_line[0].remove()

        self.previous_line = line

        self.first_plot = True

        if self.x_lim:
            self.spectre_plot.set_xlim(self.x_lim)
        else:
            self.x_lim = self.spectre_plot.get_xlim()

        if self.y_lim:
            self.spectre_plot.set_ylim(self.y_lim)
        else:
            self.y_lim = self.spectre_plot.get_ylim()

        self.grid.addWidget(self.spectre_plot, 0, 1)

    def previous_spectre(self, event):
        if self.index > 0:
            self.index -= 1
        else:
            self.index = len(self.sorted_file_names) -1
        self.plot()

    def next_spectre(self, event):
        if self.index < len( self.sorted_file_names ) + 1:
            self.index += 1
        else:
            self.index = 0
        self.plot()

    def zoom_callback(self, event):
        if self.first_plot:
            self.x_lim = event.get_xlim()
            self.y_lim = event.get_ylim()

    def create_movie(self):
        np.random.seed(19680801)
        files = []

        i = 0
        for file_name in self.sorted_file_names:
            print(file_name)
            data = h5py.File(self.output_path + file_name, 'r')['amplitude_corrected'][()]
            plt.cla()
            x = data[:, [0]]
            y = data[:, [3]]
            plt.rc('xtick', labelsize=8)
            plt.rc('ytick', labelsize=8)
            modified_julian_days = file_name.split("_")[1]
            label = Time(modified_julian_days, format='mjd' ).strftime( '%d %b %Y' )
            plt.plot(x, y, label=label)
            plt.xlim(self.x_lim)
            plt.ylim(self.y_lim)
            plt.legend(loc=1, prop={'size': 12})
            plt.xlabel("Velocity (km sec$^{-1}$)", fontsize=6)
            plt.ylabel("Flux density (Jy)", fontsize=6)
            plt.grid(True)
            fname = '_tmp%03d.png' % i
            i += 1
            print('Saving frame', fname)
            plt.savefig(fname, dpi=300, quality=100, format="png")
            files.append(fname)

        print('Making movie animation.mpg - this may take a while')
        subprocess.call("mencoder 'mf://_tmp*.png' -mf w=800:h=600:type=png:fps=10 -ovc lavc " 
                        "-lavcopts vcodec=mpeg4:mbd=2:trell:autoaspect -oac copy -o " +
                        self.source_name + "_" + get_args("line") + "_spectre_movie.mpg", shell=True)

    def center(self):
        """

        :return: None
        """
        frame_geometry = self.frameGeometry()
        centre_position = QDesktopWidget().availableGeometry().center()
        frame_geometry.moveCenter(centre_position)
        self.move(frame_geometry.topLeft())


def main():
    """

    :return: None
    """
    q_app = QApplication(sys.argv)
    application = SpectreTime(get_args("source"))
    application.show()
    application.showMaximized()
    sys.exit(q_app.exec_())


if __name__ == "__main__":
    main()
