#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
Monitoring tool
"""

import sys
import os
import argparse
import json
import numpy as np
import h5py
from PyQt5.QtWidgets import QApplication, QWidget, QDesktopWidget, \
    QGridLayout, QLabel, QLineEdit, QComboBox, QPushButton
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt
from utils.ploting_qt5 import Plot
from parsers.configparser_ import ConfigParser


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''',
                                     epilog="""Monitor.""")
    parser.add_argument("-c", "--config", help="Configuration cfg file",
                        type=str, default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 2.0')
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


class PlottingView(QWidget):
    """
    Base class for all views
    """
    def __init__(self):
        super(PlottingView, self).__init__()
        self.setWindowIcon(QIcon('viraclogo.png'))
        self.center()
        self.grid = QGridLayout()

    def add_widget(self, widget, row, colomn):
        self.grid.addWidget(widget, row, colomn)

    def center(self):
        """

        :return: None
        """
        frame_geometry = self.frameGeometry()
        centre_position = QDesktopWidget().availableGeometry().center()
        frame_geometry.moveCenter(centre_position)
        self.move(frame_geometry.topLeft())


class Monitoring(PlottingView):
    """
    Monitoring tool
    """

    def __init__(self):
        PlottingView.__init__(self)
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.setWindowTitle("Choose Source")
        choose_label = QLabel("Choose source")
        choose_line_label = QLabel("Choose line")
        self.add_widget(choose_label, 0, 0)
        self.add_widget(choose_line_label, 0, 1)
        self.source_input = QLineEdit()
        self.add_widget(self.source_input, 1, 0)
        self.source_line_input = QLineEdit()
        self.source_line_input.setText("6668")
        self.add_widget(self.source_line_input, 1, 1)
        choose_source_button = QPushButton("Ok", self)
        choose_source_button.clicked.connect(self.plot_monitoring)
        choose_source_button.setStyleSheet("background-color: green")
        self.add_widget(choose_source_button, 1, 3)
        combo_box2 = QComboBox(self)
        combo_box2.addItem("All")
        combo_box2.addItem("Not Flag")
        #comboBox2.activated[str].connect(self.setFlags)
        self.add_widget(combo_box2, 1, 2)
        self.monitoring_plot = None
        self.monitoring_view = None

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Return:
            if len(self.source_input.text()) > 1:
                self.plot_monitoring()

    def plot_monitoring(self):
        """

        :return: None
        """
        if len(self.source_input.text()) > 1:
            self.monitoring_view = MonitoringView(self.source_input.text(),
                                                  self.source_line_input.text())
            self.monitoring_view.show()


def create_label(experiment):
    return "Station is " + experiment.location.lower() + "\n" + "Date is " + \
           " ".join(experiment.Date.split("_")) + "\n" + "Iteration number " + \
           str(experiment.Iteration_number) + "\n" + "Specie " + \
           experiment.specie + "\n" + "Type " + experiment.type


class MonitoringView(PlottingView):
    """
    Monitoring View
    """

    def __init__(self, source, line):
        PlottingView.__init__(self)
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.setWindowTitle("Monitoring")
        self.source = source
        self.line = line
        self.specter_view = None
        self.new_spectre = True
        self.multiple_spectre = True
        self.specter_plots_files = set()

        result_path = get_configs("paths", "resultFilePath")
        result_file_name = source + "_" + line + ".json"

        with open(result_path + "/" + result_file_name) as result_file:
            result_data = json.load(result_file)

        experiments = [MonitoringView.Experiment(**result_data[exper]) for exper in result_data]
        experiments.sort(key=lambda e: e.modifiedJulianDays)
        labels2 = [create_label(e) for e in experiments]

        self.monitoring_plot = Plot()
        self.monitoring_plot.creatPlot(self.grid, "Time", "Flux density (Jy)",
                                       get_configs("Full_source_name", source), (1, 0), "log")

        symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        dates = [e.modifiedJulianDays for e in experiments]
        source_velocities = get_configs('velocities', self.source + "_" + self.line).split(",")
        source_velocities = [x.strip() for x in source_velocities]
        monitoring_results = [[dates]]
        self.iterations = [e.Iteration_number for e in experiments]

        for i in range(0, len(source_velocities)):
            self.monitoring_plot.plot(dates, [e.polarizationU1[i][1] for e in experiments],
                                      symbols[i] + colors[i],
                                      fontsize=8, visible=False, picker=False)
            self.monitoring_plot.plot(dates, [e.polarizationU9[i][1] for e in experiments],
                                      symbols[i] + colors[i],
                                      fontsize=8, visible=False, picker=False)
            self.monitoring_plot.plot(dates, [e.polarizationAVG[i][1] for e in experiments],
                                      symbols[i] + colors[i], fontsize=8,
                                      label="Velocity " + source_velocities[i],
                                      visible=True, picker=5)
            monitoring_results.append([e.polarizationAVG[i][1] for e in experiments])

        np.save(get_configs("paths", "monitoringFilePath") +
                self.source + "_" + str(self.line), np.transpose(monitoring_results))
        self.monitoring_plot.addCursor(labels2)
        self.monitoring_plot.addPickEvent(self.choose_spectrum)
        self.add_widget(self.monitoring_plot, 0, 0)

    class Experiment:
        """
        Experiment class
        """

        def __init__(self, **entries):
            self.__dict__.update(entries)

    def choose_spectrum(self, event):
        """

        :return: None
        """
        index = int(event.ind[0])
        iteration = self.iterations[index]
        file = [f for f in (os.listdir(get_configs("paths", "outputFilePath") + "/" + self.line))
                if str(iteration) in f and f.startswith(self.source)]
        if len(file) > 0:
            if self.new_spectre:
                output_file = get_configs("paths", "outputFilePath") \
                              + "/" + self.line + "/" + file[0]
                self.specter_plots_files.add(output_file)
                self.specter_view = SpecterView(self.specter_plots_files, self.source)
                self.specter_view.show()
                self.new_spectre = False

            else:
                output_file = get_configs("paths", "outputFilePath") \
                              + "/" + self.line + "/" + file[0]
                self.specter_plots_files.add(output_file)
                self.specter_view.set_specter_plots_files(self.specter_plots_files)

        else:
            print("No output files for iteration " + str(iteration))

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Shift:
            self.new_spectre = True

        elif event.key() == Qt.Key_Alt:
            if self.multiple_spectre:
                self.multiple_spectre = False
                print("single spectre")
            else:
                self.multiple_spectre = True
                print("multiple spectre")


class SpecterView(PlottingView):
    """
    Monitoring View
    """

    def __init__(self, spectre_files, source):
        PlottingView.__init__(self)
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.setWindowTitle("Spectre")
        self.spectre_files = list(spectre_files)
        self.source = source
        self.polarization = "AVG"
        source_name = get_configs("Full_source_name", self.source)

        self.specter_plot = Plot()
        self.specter_plot.creatPlot(self.grid, "Velocity (km sec$^{-1}$)",
                                    "Flux density (Jy)", source_name, (1, 0), "linear")
        self.specter_plot.set_tick_params(axis="x", direction="in",
                                          which="both", length=16,
                                          width=2, labelsize=12, rotation=0)
        self.specter_plot.set_tick_params(axis="y", direction="in",
                                          which="major", length=16,
                                          width=2, labelsize=12, rotation=0)
        self.specter_plot.set_tick_params(axis="y", direction="in",
                                          which="minor", length=10,
                                          width=1.5, labelsize=12, rotation=0)
        self.add_widget(self.specter_plot, 0, 0)
        self.plot()

    def plot(self):
        print("yes")
        symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
        amplitude_colon = 3
        if self.polarization == "U1":
            amplitude_colon = 1
        elif self.polarization == "U9":
            amplitude_colon = 2
        elif self.polarization == "AVG" or self.polarization == "ALL":
            amplitude_colon = 3
        for spectre_file in self.spectre_files:
            index = self.spectre_files.index(spectre_file)
            data = h5py.File(spectre_file, 'r')['amplitude_corrected'][()]
            x = data[:, 0]
            y = data[:, amplitude_colon]
            tmp = spectre_file.split("/")[-1].split(".")[0].split("_")
            plot_name = "".join([tmp[1], tmp[2], tmp[3], tmp[4]])
            self.specter_plot.plot(x, y, symbols[index], label=plot_name)

    def set_specter_plots_files(self, specter_plots_files):
        self.spectre_files.extend(specter_plots_files)
        self.spectre_files = list(set(self.spectre_files))
        self.plot()
        self.specter_plot.draw()
        self.specter_plot.canvasShow()


def main():
    """

    :return: None
    """
    q_app = QApplication(sys.argv)
    application = Monitoring()
    application.show()
    sys.exit(q_app.exec_())


if __name__ == "__main__":
    main()