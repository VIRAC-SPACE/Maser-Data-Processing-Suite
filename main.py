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
    parser.add_argument("-m", "--manual", help="Set manual log data", action='store_true')
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


def find_log_file(log_list, iteration, line):
    """

    :param log_list: list of log files
    :param iteration: iteration of observations
    :param line: frequency
    :return:
    """
    tmpl = -1
    for log in range(0, len(log_list)):
        iteration_tmp = log_list[log].split("/")[-1].split(".")[0].split("_")[-1]
        lin = log_list[log].split("/")[-1].split(".")[0].split("_")[-2]

        if lin + "_" + iteration_tmp == "f" + line + "_" + iteration:
            tmpl = log
            break
    if tmpl == -1:
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


def create_iteration_list(path, source, line):
    """

    :param line: frequency
    :param source: source
    :param path: input file path
    :return: None
    """
    iterations = [get_iteration(iteration) for iteration in os.listdir(path) if
                  source in iteration and line in iteration]
    iterations.sort(key=int, reverse=False)
    return iterations


def create_log_file_list(path, source, line):
    """

    :param line: frequency
    :param path: log file path
    :param source: source
    :return: all log files for source
    """
    return [log for log in os.listdir(path) if log.startswith(source)]


def main():
    """
    :return: None
    """
    source_name = get_args("source")
    line = get_args("line")
    data_files_path = get_configs('paths', "dataFilePath")
    result_path = get_configs('paths', "resultFilePath")
    log_path = get_configs('paths', "logPath")

    if os.path.exists(data_files_path):
        sdr_iterations = create_iteration_list(data_files_path, source_name, line)
    else:
        sdr_iterations = []

    log_path = log_path + "SDR/"
    logfile_list = create_log_file_list(log_path, source_name, line)
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

    for experiment in result:
        if experiment.split("_")[-1] in sdr_iterations and \
                experiment.split("_")[-1] not in \
                processed_iteration and result[experiment]["type"] == "SDR":
            processed_iteration.append(experiment.split("_")[-1])

        if experiment.split("_")[-1] in processed_iteration and \
                experiment.split("_")[-1] in processed_iteration and \
                result[experiment]["type"] == "SDR" and result[experiment]["flag"] == True:
            processed_iteration.remove(experiment.split("_")[-1])

    processed_iteration.sort(key=int, reverse=False)

    try:
        for i in sdr_iterations:
            if i not in processed_iteration:
                print(i, str(logfile_list [find_log_file(logfile_list, i, get_args("line"))]))
                sys.exit()
                frequency_shifting_parameter = source_name + " " + \
                                               get_args("line") + " " + i + " " + \
                                               str(logfile_list
                                                   [find_log_file(logfile_list, i,
                                                                  get_args("line"))])
                LOGGER.info("Executing python3 " + "sdr_fs.py " +
                            frequency_shifting_parameter)
                os.system("python3 " + "sdr_fs.py " +
                          frequency_shifting_parameter)
        sys.exit()
        data_files = list()
        print(sdr_path + "/" + get_args("line"))
        if os.path.exists(sdr_path + "/" + get_args("line")):
            for data in os.listdir(sdr_path + "/" + get_args("line")):
                if data.startswith(source_name) and data.endswith(".dat"):
                    data_files.append(data)

            for d in data_files:
                if d.split(".")[0].split("_")[-1] not in sdr_processed_iteration:
                    LOGGER.info("Executing python3 " +
                                "src/totalSpectrumAnalyer_qt5.py " + d + " " + get_args("line"))
                    os.system("python3 " +
                              "src/totalSpectrumAnalyer_qt5.py " + d + " " + get_args("line"))

    except IOError as e:
        print("IO Error", e)
        sys.exit(1)

    except IndexError as e:
        print("Index Error", e)
        sys.exit(1)

    except ValueError as e:
        print("Cannot crate modified Julian Days", e)

    except TypeError as e:
        print("TypeError", e)
        sys.exit(1)

    except:
        print("Unexpected error:", sys.exc_info()[0])
        sys.exit(1)


if __name__ == "__main__":
    main()
    sys.exit(0)
