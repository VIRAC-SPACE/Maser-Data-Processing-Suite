#! /usr/bin/python3
# -*- coding: utf-8 -*-
import sys
import os
import argparse

from PyQt5.QtWidgets import (QApplication)

from parsers.configparser_ import ConfigParser

from main_view import MainView


def parse_arguments():
    """
    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''automatically call sdr_fs.py 
    and total_spectrum_analyzer_qt5.py.''', epilog="""Main program.""")
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("line", help="frequency", type=int)
    parser.add_argument("-c", "--config", help="Configuration "
                                               "cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 3.0')
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


class MainGUI(QApplication):
    def __init__(self, sys_argv):
        super(MainGUI, self).__init__(sys_argv)
        source_name = get_args("source")
        line = get_args("line")
        self.main_view = MainView(source_name, line, get_args("config"))
        self.main_view.show()


def main():
    app = MainGUI(sys.argv)
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
