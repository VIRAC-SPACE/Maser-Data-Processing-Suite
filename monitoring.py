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
from matplotlib import ticker
from astropy.timeseries import LombScargle
from astropy.io import ascii
from astropy.time import Time
import matplotlib.tri as mtri
import h5py
from PyQt5.QtWidgets import QApplication, QWidget, QDesktopWidget, \
    QGridLayout, QLabel, QLineEdit, QComboBox, QPushButton, QGroupBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt
from utils.ploting_qt5 import Plot
from utils.help import find_nearest_index
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


def create_label(experiment):
    """

    :param experiment: experiment type object
    :return: label for monitoring view cursor
    """
    return "Station is " + experiment.location.lower() + "\n" + "Date is " + \
           Time(experiment.modifiedJulianDays, format='mjd').strftime('%d %b %Y') + "\n" + "Iteration number " + \
           str(experiment.Iteration_number) + "\n" \
           + "Backend type " + experiment.type


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
        """

        :param widget: widget
        :param row: row
        :param colomn: colomn
        :return: None
        """
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
        self.flag = "All"
        combo_box2 = QComboBox(self)
        combo_box2.addItem("All")
        combo_box2.addItem("Not Flag")
        combo_box2.activated[str].connect(self.set_flag)
        self.add_widget(combo_box2, 1, 2)
        self.monitoring_plot = None
        self.monitoring_view = None

    def set_flag(self, flag):
        self.flag = flag

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
                                                  self.source_line_input.text(), self.flag)
            self.monitoring_view.show()


class ChangeView(PlottingView):
    def __init__(self, component, dates, amplitude):
        PlottingView.__init__(self)
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.setWindowTitle("Change view")
        self.component = component
        self.dates = dates
        self.amplitude = amplitude
        self.changes_values = []

        self.change_plot = Plot()
        self.change_plot.creatPlot(self.grid, "Time [mjd]", "Flux density (Jy)", component, (1, 0), "linear")
        self.change_plot.graph.scatter(self.dates, self.amplitude, color="k", s=0.01)

        changes = self.compute_changes()
        changes_size = np.abs(changes)
        change_color = []

        for change in changes:
            if change >= 0:
                change_color.append("g")
            else:
                change_color.append("r")

        self.change_plot.graph.scatter(self.dates, self.amplitude, marker="s", color=change_color, s=changes_size)
        self.change_plot.addCursor(np.round(changes,  decimals=3))
        self.change_plot.addClickEvent(self.compute)

        self.add_widget(self.change_plot, 0, 0)

    def compute_changes(self):
        changes = [0]
        for i in range(1, len(self.amplitude)):
            changes.append(((self.amplitude[i] - self.amplitude[i-1])/self.amplitude[i])*100)
        return changes

    def compute(self, event):
            if len(self.changes_values) == 0 or len(self.changes_values) == 1:
                self.changes_values.append(event.ydata)
                if len(self.changes_values) == 2:
                    print("changes",((self.changes_values[0] - self.changes_values[1]) / self.changes_values[0]) * 100)
                    self.changes_values[0] = self.changes_values[1]
                    self.changes_values[1] = event.ydata

            elif len(self.changes_values) == 2:
                print("changes", ((self.changes_values[0] - self.changes_values[1])/self.changes_values[0])*100)
                self.changes_values[0] = self.changes_values[1]
                self.changes_values[1] = event.ydata


class MonitoringView(PlottingView):
    """
    Monitoring View
    """

    def __init__(self, source, line, flag):
        PlottingView.__init__(self)
        self.min_mjd = None
        self.max_mjd = None
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.setWindowTitle("Monitoring")
        self.source = source
        self.line = line
        self.flag = flag
        self.polarization = "polarization AVG"
        self.specter_view = None
        self.period_view = None
        self.maps_view = None
        self.component_input = None
        self.change_view = None
        self.add_widget(self.create_control_group(), 0, 1)
        self.new_spectre = True
        self.multiple_spectre = True
        self.spectrum_set = set()
        self.specter_plots_files = set()
        self.lines = []
        self.flagged_points = []
        self.flags = []
        self.un_flags = []

        result_path = get_configs("paths", "resultFilePath")
        result_file_name = self.source + "_" + self.line + ".json"

        with open(result_path + "/" + result_file_name) as result_file:
            result_data = json.load(result_file)

        self.experiments = [MonitoringView.Experiment(**result_data[experiment])
                            for experiment in result_data]
        if self.flag == "Not Flag":
            self.experiments = [e for e in self.experiments if not e.flag]
        self.experiments.sort(key=lambda e: float(e.modifiedJulianDays))
        labels2 = [create_label(e) for e in self.experiments]

        self.monitoring_plot = Plot()
        self.monitoring_plot.creatPlot(self.grid, "Time", "Flux density (Jy)",
                                       get_configs("Full_source_name", source), (1, 0), "log")

        symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        self.dates = [float(e.modifiedJulianDays) for e in self.experiments]
        self.source_velocities = get_configs('velocities', self.source + "_" + self.line).split(",")
        self.source_velocities = [x.strip() for x in self.source_velocities]
        monitoring_results = [np.array([self.dates])]
        self.iterations = [e.Iteration_number for e in self.experiments]

        self.line_dict = {"left": list(),
                          "right": list(),
                          "avg": list()}

        for e in self.experiments:
            if len(e.polarizationAVG) != len(self.source_velocities):
                print(e.Iteration_number)

        print(self.experiments[0].areas)
        for i in range(0, len(self.source_velocities)):
            l1 = self.monitoring_plot.plot(self.dates,
                                           [e.polarizationU1[i][1] for e in self.experiments],
                                           symbols[i] + colors[i],
                                           fontsize=8, visible=False, picker=False)
            l2 = self.monitoring_plot.plot(self.dates,
                                           [e.polarizationU9[i][1] for e in self.experiments],
                                           symbols[i] + colors[i],
                                           fontsize=8, visible=False, picker=False)
            l3 = self.monitoring_plot.plot(self.dates,
                                           [e.polarizationAVG[i][1] for e in self.experiments],
                                           symbols[i] + colors[i], fontsize=8,
                                           label="Velocity " + self.source_velocities[i],
                                           visible=True, picker=5)
            monitoring_results.append(np.array([e.polarizationAVG[i][1] for e in self.experiments]))

            self.lines.append(l1)
            self.lines.append(l2)
            self.lines.append(l3)
            self.line_dict["left"].append(l1)
            self.line_dict["right"].append(l2)
            self.line_dict["avg"].append(l3)

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
            self.flag = None
            self.modifiedJulianDays = None
            self.Iteration_number = None
            self.polarizationAVG = None
            self.polarizationU9 = None
            self.polarizationU1 = None
            self.location = None
            self.Date = None
            self.specie = None
            self.time = None
            self.type = None
            self.__dict__.update(entries)

    def choose_spectrum(self, event):
        """

        :return: None
        """
        try:
            ind = event.ind[0]
        except AttributeError as error:
            print("Attribute Error", error, sys.exc_info()[0])
        else:
            this_line = event.artist
            xdata = this_line.get_xdata()
            ydata = this_line.get_ydata()
            result_path = get_configs("paths", "resultFilePath")
            result_file_name = self.source + "_" + self.line + ".json"
            iteration = self.iterations[ind]
            mjd = [e.modifiedJulianDays for e in self.experiments][ind]
            station = [e.location for e in self.experiments][ind]

            with open(result_path + result_file_name) as result_data:
                results = json.load(result_data)

            if event.mouseevent.button == 1:
                output_file = get_configs("paths", "outputFilePath") + self.line + "/" + self.source + "/" + \
                              self.source + "_" + str(mjd) + "_" + \
                              str(station) + "_" + str(iteration) + ".h5"
                if os.path.isfile(output_file):
                    if self.new_spectre:
                        self.specter_plots_files.clear()
                        self.specter_plots_files.add(output_file)
                        self.specter_view = SpecterView(self.specter_plots_files,
                                                        self.source, self.polarization)
                        self.spectrum_set.add(self.specter_view)
                        self.specter_view.show()
                        self.new_spectre = False

                    else:
                        self.specter_plots_files.add(output_file)
                        self.specter_view.set_specter_plots_files(self.specter_plots_files)

                else:
                    print("No output files for iteration " + str(iteration))

            elif event.mouseevent.button == 2:
                if len(self.flagged_points) > 0:
                    if (xdata[ind], ydata[ind]) not in self.un_flags:
                        selected_point = (xdata[ind], ydata[ind])
                        flagged_points_data = list(set([(point.get_xdata()[0],
                                                         point.get_ydata()[0]) for
                                                        point in self.flagged_points]))
                        unflag_index = find_nearest_index(flagged_points_data, selected_point)
                        if unflag_index != 0:
                            unflag_index -= 1

                        self.flagged_points[unflag_index].set_visible(False)
                        self.flagged_points[unflag_index].remove()
                        self.flagged_points.pop(unflag_index)
                        self.flags.pop(unflag_index)
                        self.monitoring_plot.canvasShow()

                        for experiment in results:
                            if experiment.endswith("_" + str(iteration)) and experiment.split("_")[1] == str(mjd):
                                results[experiment]["flag"] = False
                                print(experiment, "is un flag")

                        with open(get_configs("paths", "resultFilePath") +
                                  result_file_name, "w") as result_data:
                            result_data.write(json.dumps(results, indent=2))

                        self.un_flags.append((xdata[ind], ydata[ind]))

            elif event.mouseevent.button == 3:
                if (xdata[ind], ydata[ind]) not in self.flags:
                    flagged_point = self.monitoring_plot.plot(xdata[ind],
                                                              ydata[ind], "rx", markersize=10)[0]
                    self.flagged_points.append(flagged_point)
                    self.monitoring_plot.canvasShow()

                    if len(self.un_flags) > 0:
                        index_tmp = find_nearest_index(self.un_flags, (xdata[ind], ydata[ind]))
                        if index_tmp > 0:
                            index_tmp -= 1
                        self.un_flags.pop(index_tmp)

                    for experiment in results:
                        if experiment.endswith("_" + str(iteration)) and experiment.split("_")[1] == str(mjd):
                            results[experiment]["flag"] = True
                            print(experiment, "is flag")

                    with open(get_configs("paths", "resultFilePath") +
                              result_file_name, "w") as result_data:
                        result_data.write(json.dumps(results, indent=2))
                    self.flags.append((xdata[ind], ydata[ind]))

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

    def set_polarization(self, polarization):
        self.polarization = polarization
        self.change_visible_lines(polarization)

    def create_control_group(self):
        group_box = QGroupBox("")
        group_box.setFixedWidth(120)
        group_box.setFixedHeight(150)
        control_grid = QGridLayout()
        polarization_combo_box = QComboBox(self)
        polarization_combo_box.addItem("polarization AVG")
        polarization_combo_box.addItem("polarization left")
        polarization_combo_box.addItem("polarization right")
        polarization_combo_box.addItem("ALL")
        polarization_combo_box.activated[str].connect(self.set_polarization)
        control_grid.addWidget(polarization_combo_box, 5, 0)
        plot_periods_button = QPushButton('Plot periods', self)
        plot_periods_button.clicked.connect(self.create_period_view)
        control_grid.addWidget(plot_periods_button, 3, 0)
        self.max_mjd = QLineEdit()
        self.max_mjd.setFixedWidth(100)
        control_grid.addWidget(self.max_mjd, 0, 0)
        self.min_mjd = QLineEdit()
        self.min_mjd.setFixedWidth(100)
        control_grid.addWidget(self.min_mjd, 1, 0)
        self.component_input = QLineEdit()
        self.component_input.setFixedWidth(100)
        control_grid.addWidget(self.component_input, 2, 0)
        plot_maps_button = QPushButton('Plot maps', self)
        plot_maps_button.clicked.connect(self.create_map_view)
        control_grid.addWidget(plot_maps_button, 4, 0)

        plot_changes_button = QPushButton('Plot changes', self)
        plot_changes_button.clicked.connect(self.create_change_view)
        control_grid.addWidget(plot_changes_button, 5, 0)

        group_box.setLayout(control_grid)
        return group_box

    def change_visible_lines(self, polarization):
        if polarization == "ALL":
            for line in self.lines:
                line[0].set_picker(5)
                line[0].set_visible(True)

        elif polarization == "polarization AVG":
            for line in self.line_dict["avg"]:
                line[0].set_picker(5)
                line[0].set_visible(True)

            for line in self.line_dict["left"]:
                line[0].set_picker(False)
                line[0].set_visible(False)

            for line in self.line_dict["right"]:
                line[0].set_picker(False)
                line[0].set_visible(False)

        elif polarization == "polarization left":
            for line in self.line_dict["left"]:
                line[0].set_picker(5)
                line[0].set_visible(True)

            for line in self.line_dict["avg"]:
                line[0].set_picker(False)
                line[0].set_visible(False)

            for line in self.line_dict["right"]:
                line[0].set_picker(False)
                line[0].set_visible(False)

        elif polarization == "polarization right":
            for line in self.line_dict["right"]:
                line[0].set_picker(5)
                line[0].set_visible(True)

            for line in self.line_dict["left"]:
                line[0].set_picker(False)
                line[0].set_visible(False)

            for line in self.line_dict["avg"]:
                line[0].set_picker(False)
                line[0].set_visible(False)

        self.monitoring_plot.canvasShow()

    def create_period_view(self):
        component = self.component_input.text()
        max_mjd = self.max_mjd.text()
        min_mjd = self.min_mjd.text()
        if component in self.source_velocities:
            component_index = self.source_velocities.index(component)
            amplitude = [e.polarizationAVG[component_index][1] for e in self.experiments]
            symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
            colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
            plot_simbol = symbols[component_index] + colors[component_index]
            self.period_view = PeriodView(self.dates, amplitude, plot_simbol, component, max_mjd, min_mjd)
            self.period_view.show()
        else:
            print("wrong velocity selected")

    def create_map_view(self):
        self.maps_view = MapsView(self.dates, self.source, self.line)
        self.maps_view.show()

    def create_change_view(self):
        component = self.component_input.text()
        if component in self.source_velocities:
            component_index = self.source_velocities.index(component)
            amplitude = [e.polarizationAVG[component_index][1] for e in self.experiments]
            component = self.component_input.text()
            self.change_view = ChangeView(component, self.dates, amplitude)
            self.change_view.show()


class SpecterView(PlottingView):
    """
    Monitoring View
    """

    def __init__(self, spectre_files, source, polarization):
        PlottingView.__init__(self)
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.setWindowTitle("Spectre")
        self.spectre_files = list(set(spectre_files))
        self.source = source
        self.polarization = polarization
        self.plot_set = set()
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
        symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
        amplitude_colon = 3
        if self.polarization == "polarization left":
            amplitude_colon = 1
        elif self.polarization == "polarization right":
            amplitude_colon = 2
        elif self.polarization == "polarization AVG" or self.polarization == "ALL":
            amplitude_colon = 3
        for spectre_file in self.spectre_files:
            index = self.spectre_files.index(spectre_file)
            data_file = h5py.File(spectre_file, 'r')
            if "amplitude_corrected_not_smooht" in data_file:
                data = data_file['amplitude_corrected_not_smooht'][()]
                x = data[:, 0]
                y = data[:, amplitude_colon]
                plot_name = ".".join([spectre_file.split("/")[-1].split(".")[0],
                                      spectre_file.split("/")[-1].split(".")[1]])
                if plot_name not in self.plot_set:
                    self.specter_plot.plot(x, y, symbols[index], label=plot_name)
                    self.plot_set.add(plot_name)
            else:
                print("Output " + spectre_file + " file has no amplitude_corrected_not_smooht colomm")

    def set_specter_plots_files(self, specter_plots_files):
        """

        :param specter_plots_files: total_spectrum_analyser_qt5.py output files
        :return: None
        """
        self.spectre_files.extend(specter_plots_files)
        self.spectre_files = list(set(self.spectre_files))
        self.plot()
        self.specter_plot.draw()
        self.specter_plot.canvasShow()


class PeriodView(PlottingView):
    """
    Period View
    """
    def __init__(self, time, amplitude, plot_symbol, velocity_name, max_mjd, min_mjd):
        PlottingView.__init__(self)
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.setWindowTitle("Periods in days ")
        self.time = time
        self.amplitude = amplitude
        self.plot_symbol = plot_symbol
        self.velocity_name = velocity_name
        self.max_mjd = max_mjd
        self.min_mjd = min_mjd

        if self.max_mjd != "" or self.min_mjd != "":
            self.max_mjd_index = find_nearest_index(self.time, float(self.max_mjd))
            self.min_mjd_index = find_nearest_index(self.time, float(self.min_mjd))
            self.time = self.time[self.min_mjd_index:self.max_mjd_index]
            self.amplitude = self.amplitude[self.min_mjd_index:self.max_mjd_index]

        error = np.array(self.amplitude) * 0.1
        period_range = [10, 100]
        ls_range = np.linspace(period_range[0], period_range[1], 10000)
        freq = 1. / ls_range
        time = np.array(self.time) - self.time[0]

        ls = LombScargle(time, self.amplitude, error, fit_mean=True, normalization='model', nterms=1)
        power = ls.power(freq)
        times = 1. / freq

        self.period_plot = Plot()
        self.period_plot.creatPlot(self.grid, "Period (days)", "Power", None, (1, 0), "linear")
        self.add_widget(self.period_plot, 0, 0)
        self.period_plot.plot(times, power, self.plot_symbol,
                              label="polarization AVG " + "Velocity " + self.velocity_name, rasterized=True)


class MapsView(PlottingView):
    """
    Maps View
    """
    def __init__(self, mjd, source, line):
        PlottingView.__init__(self)
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.setWindowTitle(" ")
        self.mjd = mjd
        self.source = source
        self.line = line
        self.output_files = os.listdir(get_configs("paths", "outputFilePath") + self.line + "/" + self.source)
        self.output_files = [of for of in self.output_files if of.startswith(self.source)]

        days = self.mjd[-1] - self.mjd[0]
        sources_vrange = ascii.read('DB_vrange.csv')
        source_vrange_index = sources_vrange['name'].tolist().index(self.source)
        vmin = dict(sources_vrange)["vmin"][source_vrange_index]
        vmax = dict(sources_vrange)["vmax"][source_vrange_index]

        if vmin is None or vmax is None:
            vmin, vmax, velocity, observed_flux, observed_time = self.get_max_min_velocity(vmin, vmax)

        else:
            _, _, velocity, observed_flux, observed_time = self.get_max_min_velocity(vmin, vmax)

        observed_flux = list(np.array(observed_flux).clip(min=1.5))
        triang = mtri.Triangulation(velocity, observed_time)

        self.map_plot = Plot()
        self.map_plot.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', "MJD", None, (1, 0), "log")
        lvls = np.linspace(int(np.min(observed_flux)), int(np.max(observed_flux)), 1000)
        cs = self.map_plot.graph.tricontourf(triang, observed_flux, levels=lvls,
                                             antialiased=False, locator=ticker.LogLocator, cmap="jet")

        cs.set_clim(vmin=1.5)
        cbar = self.map_plot.colorbar(cs, spacing="proportional", label=r'$Flux~(\mathrm{Jy})$', extendrect=False)
        cbar.locator = ticker.LogLocator()
        self.add_widget(self.map_plot, 0, 0)

        self.map_plot.save_fig("cepa_maps.pdf", format="pdf", dpi=150)

    def get_max_min_velocity(self, vmin, vmax):
        max_velocitys = []
        min_velocitys = []
        velocity_ = []
        observed_flux = []
        observed_time = []

        for file in self.output_files:
            data_file = h5py.File(get_configs("paths", "outputFilePath") +
                                  self.line + "/" +
                                  self.source + "/" + file, "r")
            if "amplitude_corrected_not_smooht" in data_file:
                data = data_file["amplitude_corrected_not_smooht"][()]
                velocity = data[:, 0]
                amplitude = data[:, 3]
                mjd = file.split("_")[1]

                max_velocitys.append(max(velocity))
                min_velocitys.append(min(velocity))

                for i in range(0, len(velocity)):
                    if vmin or vmax:
                        if vmin <= velocity[i] <= vmax:
                            velocity_.append(velocity[i])
                            observed_flux.append(amplitude[i])
                            observed_time.append(float(mjd))
                    else:
                        velocity_.append(velocity[i])
                        observed_flux.append(amplitude[i])
                        observed_time.append(float(mjd))

        return max(max_velocitys), min(min_velocitys), velocity_, observed_flux, observed_time


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