#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
run scripts sdr_fs.py and total_spectrum_analyzer_qt5.py for all unprocessed experiments
"""
import os
import sys
import argparse
import json
import logging
import coloredlogs
import h5py
from parsers.configparser_ import ConfigParser

coloredlogs.install(level='PRODUCTION')
LOGGER = logging.getLogger('Main')


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


def find_log_file(log_list, iteration, station):
    """

    :param log_list: list of log files
    :param iteration: iteration of observations
    :param station: station of observation
    :return:
    """
    tmpl = ""
    for log in log_list:
        if "_" + str(iteration) in log and station in log:
            tmpl = log
            break

    else:
        tmpl = log_list[-1]

    if tmpl == "":
        LOGGER.warning("Warning " + "log for iteration " +
                       iteration + " do not exist log file " +
                       log_list[-1] + " will be used instead!")
    return tmpl


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
    iterations = [get_iteration(iteration) for iteration in os.listdir(path) if
                  source in iteration and line in iteration and os.path.isdir(path + iteration)]
    iterations.sort(key=int, reverse=False)
    return iterations


def create_station_list(path, source, line):
    """

    :param line: frequency
    :param source: source
    :param path: input file path
    :return: stations list
    """
    iterations = [file for file in os.listdir(path) if source in file and line in file and os.path.isdir(path + file)]
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


def main():
    """
    :return: None
    """
    source_name = get_args("source")
    line = get_args("line")
    data_files_path = get_configs('paths', "dataFilePath")
    result_path = get_configs('paths', "resultFilePath")
    log_path = get_configs('paths', "logPath")
    output_path = get_configs('paths', "outputFilePath")

    if os.path.exists(data_files_path):
        sdr_iterations = create_iteration_list(data_files_path, source_name, line)
    else:
        sdr_iterations = []

    log_path = log_path + "SDR/"
    logfile_list = create_log_file_list(log_path, source_name, line)
    station_list = create_station_list(data_files_path, source_name, line)
    result_file_name = result_path + source_name + "_" + get_args("line") + ".json"

    if os.path.isfile(result_file_name):
        pass
    else:
        os.system("touch " + result_file_name)
        result_file = open(result_file_name, "w")
        result_file.write("{ \n" + "\n}")
        result_file.close()

    with open(result_file_name, "r") as result_data:
        result = json.load(result_data)

    processed_iteration = list()
    processed_iteration2 = list()

    for experiment in result:
        if experiment.split("_")[-1] in sdr_iterations and \
                experiment.split("_")[-1] not in \
                processed_iteration and result[experiment]["type"] == "SDR":
            processed_iteration.append(experiment.split("_")[-1])

        if experiment.split("_")[-1] in processed_iteration and \
                result[experiment]["type"] == "SDR" and result[experiment]["flag"]:
            processed_iteration.remove(experiment.split("_")[-1])

        if experiment.split("_")[-1] not in processed_iteration2 and result[experiment]["type"] == "SDR":
            processed_iteration2.append(experiment.split("_")[-1])

    processed_iteration.sort(key=int, reverse=False)
    processed_iteration2.sort(key=int, reverse=False)

    station_index = 0
    for iteration in sdr_iterations:
        if iteration not in processed_iteration:
            station = station_list[station_index]
            log_file = find_log_file(logfile_list, iteration, station)
            sdr_fs_parameter = source_name + " " + line + " " + iteration + " " + log_file
            LOGGER.info("Executing python3 " + "sdr_fs.py " + sdr_fs_parameter)
            os.system("python3 " + "sdr_fs.py " + sdr_fs_parameter)
            station_index += 1

    output_files = os.listdir(output_path + "/" + line + "/" + source_name)
    for output_file in output_files:
        if output_file.split("_")[-1].split(".")[0] not in processed_iteration2:
            if output_file.startswith(source_name):
                with h5py.File(get_configs("paths", "outputFilePath") + get_args("line") +
                                "/" + get_args("source") + "/" + output_file, "r") as input_data_file:
                    input_file_keys = list( input_data_file.keys())
                if "amplitude" in input_file_keys:
                    LOGGER.info("Executing python3 " +
                                 "total_spectrum_analyzer_qt5.py " + output_file + " " + line)
                    os.system("python3 " +
                               "total_spectrum_analyzer_qt5.py " + output_file + " " + line)


if __name__ == "__main__":
    main()
    sys.exit( 0 )
