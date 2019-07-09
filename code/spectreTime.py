#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import matplotlib
import datetime

from PyQt5.QtWidgets import QWidget, QGridLayout, QApplication, QDesktopWidget, QToolButton, QLabel
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt

from ploting_qt5 import  Plot
from parsers._configparser import ConfigParser
from help import *
from result import  Result
from monitor.months import Months
from monitor.monitoringViewHelper import MonitoringViewHelper

matplotlib.use('Qt5Agg')


def parseArguments():
    parser = argparse.ArgumentParser(description='''Monitoring spectre change in time. ''', epilog="""Monitor Spectre.""")
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    return args


def getArgs(key):
    return str(parseArguments().__dict__[key])


def getConfigs(key, value):
    configFilePath = getArgs("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)


class SpectreTime(QWidget):
    def __init__(self, source):
        super().__init__()
        self.source = source
        self.source_name = getConfigs("Full_source_name", self.source)
        self.output_path = getConfigs("paths", "notSmoohtFilePath") + "/" + self.source + "/"
        output_files = [file for file in os.listdir(self.output_path) if file.startswith(self.source + "_")]
        self.index = 0
        self.months = Months()
        dates = [self.__create_date(file_name) for file_name in output_files]
        dates = sorted(dates)

        self.sorted_file_names = []

        for f in range(0, len(dates)):
            for file_name in output_files:
                if self.__create_date(file_name) == dates[f]:
                    self.sorted_file_names.append(file_name)


        self.setWindowIcon(QIcon('viraclogo.png'))
        self.center()
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.setWindowTitle(self.source_name)
        self.InformationLb = QLabel("")
        self.__addWidget(self.InformationLb, 0, 2)
        self.previous = QToolButton(arrowType=Qt.LeftArrow)
        self.next = QToolButton(arrowType=Qt.RightArrow)
        self.previous.clicked.connect(self.previous_spectre)
        self.next.clicked.connect(self.next_spectre)
        self.__addWidget(self.previous, 1, 0)
        self.__addWidget(self.next, 1, 2)
        self.plot()

    def __addWidget(self, widget, row, colomn):
        self.grid.addWidget(widget, row, colomn)


    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def __create_date(self, file_name):
        tmpDate = file_name.split("/")[-1].split("_")
        tmpDate[-4] = self.months.getMonthNumber([tmpDate[-4]][0])
        date = datetime.datetime.strptime(" ".join(tmpDate[1:-2]), "%H %M %S %d %m %Y")
        return date


    def next_spectre(self, event):
        if self.index < len(self.sorted_file_names) +1:
            self.index += 1
        else:
            self.index = 0
        self.plot()


    def previous_spectre(self, event):
        if self.index > 0:
            self.index -= 1
        else:
            self.index = len(self.sorted_file_names) -1
        self.plot()


    def plot(self):
        file_name = self.output_path  + self.sorted_file_names [self.index]
        data = np.fromfile(file_name, dtype="float64", count=-1, sep=" ").reshape((file_len(file_name), 4))

        date = self.__create_date(file_name)
        iteration = file_name.split("/")[-1].split("_")[-1].split(".")[0]
        self.InformationLb.setText("Date " + str(date) + "\n" + "Iteration " + str(iteration))

        x = data[:, [0]]
        y = data[:, [3]]

        plot = Plot()
        plot.creatPlot(self.grid, "Velocity (km sec$^{-1}$)", "Flux density (Jy)", getConfigs("Full_source_name", self.source), (1, 1), "linear")
        plot.plot(x, y, "-")
        self.__addWidget(plot, 0,1)


class Main():

    def __parse_arguments(self):
        self.source = getArgs("source")

    def __create_app(self):
        qApp = QApplication(sys.argv)
        aw = SpectreTime(self.source)
        aw.show()
        sys.exit(qApp.exec_())
        sys.exit(0)

    @classmethod
    def run(cls):
        Main.__parse_arguments(cls)
        Main.__create_app(cls)


def main():
    Main().run()
    sys.exit(0)


if __name__ == "__main__":
    main()
