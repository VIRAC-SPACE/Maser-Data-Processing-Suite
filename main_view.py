import sys
import os
import json
import logging
from configparser import NoSectionError

import coloredlogs
import h5py
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtCore import QRunnable, QThreadPool, QCoreApplication

from main_gui_view import Ui_MainGui
from parsers.configparser_ import ConfigParser

coloredlogs.install(level='PRODUCTION')
LOGGER = logging.getLogger('Main')


def get_configs(section, key, config_file):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :param config_file: configuration file path
    :return: configuration file section key value
    """
    config_file_path = config_file
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def get_iteration(dir_name):
    """

    :param dir_name:
    :return: iteration number
    """
    return dir_name.split("_")[-1]


def get_station(dir_name):
    """

    :param dir_name:
    :return: station
    """
    return dir_name.split("_")[-2]


def create_iteration_list(path, source, line):
    """

    :param line: frequency
    :param source: source
    :param path: input file path
    :return: iterations list
    """
    stations = list(set(create_station_list(path, source, line)))
    iterations_for_source_and_line = [file for file in os.listdir(path)
                                      if source + "_" in file and line in file and os.path.isdir(path + file)]
    iterations_for_station = {station: [] for station in stations}
    for iteration in iterations_for_source_and_line:
        iterations_for_station[get_station(iteration)].append(get_iteration(iteration))

    for station in stations:
        iterations_for_station[station].sort(key=int, reverse=False)

    return iterations_for_station


def create_station_list(path, source, line):
    """

    :param line: frequency
    :param source: source
    :param path: input file path
    :return: stations list
    """
    iterations = [file for file in os.listdir(path) if
                  source + "_" in file and line in file and os.path.isdir(path + file)]
    iterations.sort(key=get_iteration, reverse=False)
    stations = [get_station(iteration) for iteration in iterations]
    return stations


def create_log_file_list(path, source, line):
    """

    :param line: frequency
    :param path: log file path
    :param source: source
    :return: all log files for source
    """
    return [log for log in os.listdir(path) if log.startswith(source + "_") and line in log]


class RunSDR(QRunnable):
    def __init__(self, source_name, line, iteration, log_file, progress, update):
        super().__init__()
        self.source_name = source_name
        self.line = line
        self.iteration = iteration
        self.log_file = log_file
        self.progress = progress
        self.update = update

    def run(self):
        sdr_fs_parameter = self.source_name + " " + self.line + " " + self.iteration + " " + self.log_file
        LOGGER.info("Executing python3 " + "sdr_fs.py " + sdr_fs_parameter)
        try:
            if not os.path.exists(get_configs("paths", "logPath", "config/config.cfg") + "SDR/" + self.log_file):
                LOGGER.warning("Warning log file " + self.log_file + " do not exist")
        except NoSectionError:
            pass
        os.system("python3 " + "sdr_fs.py " + self.source_name + " " + str(self.line) + " " +
                  str(self.iteration) + " " + self.log_file)
        self.progress.setValue(self.update)


class RunTotal(QRunnable):
    def __init__(self, output_file, line, source_name, progress, update):
        super().__init__()
        self.output_file = output_file
        self.line = line
        self.source_name = source_name
        self.progress = progress
        self.update = update

    def run(self):
        total_parameter = self.output_file + " " + self.line
        LOGGER.info("Executing python3 " + "total_spectrum_analyzer_qt5.py " + total_parameter)
        os.system("python3 " + "total_spectrum_analyzer_qt5.py " + self.output_file + " " + str(self.line))
        self.progress.setValue(self.update)


class MainView(QMainWindow):
    def __init__(self, source_name, line, config, *args, **kwargs):
        super(MainView, self).__init__(*args, **kwargs)
        self._ui = Ui_MainGui()
        self._ui.setupUi(self)
        self.setWindowTitle('MainGUI')
        self.setWindowIcon(QIcon('viraclogo.png'))
        self.source_name = source_name
        self.line = line
        self.config = config
        self._ui.source_label.setText("Source: " + self.source_name)
        self._ui.line_label.setText("Line: " + self.line)
        self._ui.sdr_ProgressBar.setValue(0)
        self._ui.toral_ProgressBar.setValue(0)
        self._ui.run_button.clicked.connect(self.run)
        self._ui.quit_button.clicked.connect(self.quit)

        self.run_sdr = None
        self.sdr_pool = None
        self.total_pool = None
        self.run_total = None

    def run(self):
        data_files_path = get_configs('paths', "dataFilePath", self.config)
        result_path = get_configs('paths', "resultFilePath", self.config)
        output_path = get_configs('paths', "outputFilePath", self.config)

        if os.path.exists(data_files_path):
            sdr_iterations = create_iteration_list(data_files_path, self.source_name, self.line)
        else:
            sdr_iterations = []

        result_file_name = result_path + self.source_name + "_" + self.line + ".json"

        if os.path.isfile(result_file_name):
            pass
        else:
            os.system("touch " + result_file_name)
            result_file = open(result_file_name, "w")
            result_file.write("{ \n" + "\n}")
            result_file.close()

        with open(result_file_name, "r") as result_data:
            result = json.load(result_data)
        stations = list(set(create_station_list(data_files_path, self.source_name, self.line)))
        processed_iteration = {station: [] for station in stations}
        processed_iteration2 = {station: [] for station in stations}

        for experiment in result:
            if get_station(experiment) == "IRBENE16":
                station = "ib"
            else:
                station = "ir"

            iteration_in_result = experiment.split("_")[-1]
            if iteration_in_result in sdr_iterations[station] and \
                    iteration_in_result not in \
                    processed_iteration[station] and result[experiment]["type"] == "SDR":
                processed_iteration[station].append(get_iteration(experiment))

            if iteration_in_result in processed_iteration[station] and \
                    result[experiment]["type"] == "SDR" and result[experiment]["flag"]:
                processed_iteration[station].remove(iteration_in_result)

            if iteration_in_result not in processed_iteration2[station] and result[experiment]["type"] == "SDR":
                processed_iteration2[station].append(iteration_in_result)
        for station in stations:
            processed_iteration[station].sort(key=int, reverse=False)
            processed_iteration2[station].sort(key=int, reverse=False)

        run_sdr_iterations = []
        for station in stations:
            for iteration in sdr_iterations[station]:
                if iteration not in processed_iteration[station]:
                    run_sdr_iterations.append(iteration)

        if len(run_sdr_iterations) == 0:
            print("No SDR iterations")
        else:
            self.sdr_pool = QThreadPool.globalInstance()
            self.sdr_pool.setMaxThreadCount(1)
            sdr_count = 0
            for station in stations:
                for iteration in sdr_iterations[station]:
                    if iteration not in processed_iteration[station]:
                        log_file = self.source_name + "_" + "f" + self.line + "_" + station + "_" + iteration + ".log"
                        sdr_count += 1
                        self.run_sdr = RunSDR(self.source_name, self.line, iteration, log_file, self._ui.sdr_ProgressBar,
                                              sdr_count/len(run_sdr_iterations) * 100)
                        self.sdr_pool.start(self.run_sdr)

        if self.sdr_pool:
            self.sdr_pool.waitForDone()

        output_files = os.listdir(output_path + "/" + self.line + "/" + self.source_name)
        run_total_iterations = []
        for output_file in output_files:
            output_file_station = output_file.split("_")[-2].split(".")[0]
            output_file_iteration = output_file.split("_")[-1].split(".")[0]
            if output_file_station == "IRBENE16":
                station = "ib"
            else:
                station = "ir"

            if output_file_iteration not in processed_iteration2[station]:
                if output_file.startswith(self.source_name):
                    with h5py.File(get_configs("paths", "outputFilePath", "config/config.cfg") +
                                   str(self.line) + "/" + self.source_name + "/" + output_file, "r") as input_data_file:
                        input_file_keys = list(input_data_file.keys())
                        input_data_file.close()
                        if "amplitude" in input_file_keys:
                            run_total_iterations.append(output_file)

        if len(run_total_iterations) == 0:
            print("No TOTAL iterations")
        else:
            self.total_pool = QThreadPool.globalInstance()
            self.total_pool.setMaxThreadCount(1)
            total_count = 0
            for output_file in run_total_iterations:
                total_count += 1
                self.run_total = RunTotal(output_file, self.line, self.source_name, self._ui.toral_ProgressBar,
                                          total_count/len(run_total_iterations) * 100)
                self.total_pool.start(self.run_total)

    def quit(self):
        print("quit")
        try:
            if self.run_sdr and self.sdr_pool:
                self.sdr_pool.cancel(self.run_sdr)
                self.run_sdr.autoDelete()
                del self.run_sdr
                self.sdr_pool.clear()
                del self.sdr_pool
            if self.run_sdr:
                self.run_sdr.autoDelete()
                del self.run_sdr
            if self.sdr_pool:
                self.sdr_pool.clear()
                del self.sdr_pool
            if self.run_total:
                self.run_total.autoDelete()
                del self.run_total
            if self.total_pool:
                self.total_pool.clear()
                del self.total_pool
        except AttributeError:
            pass
        except RuntimeError:
            pass

        self.close()
        self.destroy()
        del self
        QCoreApplication.instance().quit()
        sys.exit(1)
