import sys
import os
import argparse
from functools import reduce
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser
from utils.help import file_len, correct_numpy_read_data, convert_datetime_object_to_mjd


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''c. ''')
    parser.add_argument("line", help="line", type=int, default=6668)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="../config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 2.0')
    args = parser.parse_args()
    return args


def get_args(key):
    """

    :param key: argument key
    :return: to script passed argument value
    """
    return str(parse_arguments().__dict__[key])


def get_configs_items(section):
    """

    :return: None
    """
    config_file_path = get_args("config")
    config = ConfigParser(config_file_path)
    return config.get_items(section)


def get_configs(section, key):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :return: configuration file section key value
    """
    config_file_path = get_args("config")
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def main():
    old_monitoring_file_path = get_configs("paths", "oldMonitoringFilePath")
    new_monitoring_file_path = get_configs("paths", "monitoringFilePath")

    all_sources = list(get_configs_items("sources").keys())

    function_indexes = []
    variability_indexes = []
    mean_of_y = []
    for source in all_sources:
        new_monitoring_file = new_monitoring_file_path + "/" + source + "_" + get_args("line") + ".npy"
        old_monitoring_file = old_monitoring_file_path + "/" + source + ".dat"

        if os.path.isfile(old_monitoring_file) and os.path.isfile(new_monitoring_file):
            component_count = len(get_configs("velocities",
                                              source + "_" +
                                              get_args("line")).replace(" ", "").split(","))
            components = [i for i in range(1, component_count + 1)]
            new_data = np.load(new_monitoring_file, allow_pickle=True)
            old_data = np.loadtxt(old_monitoring_file, dtype=str).reshape(
                (file_len(old_monitoring_file), component_count + 1))
            old_x = correct_numpy_read_data(old_data[:, [0]])
            old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
            old_data[:, [0]] = old_x[0]
            old_data_tmp = []
            for tmp in range(0, len(new_data)):
                old_data_tmp.append([])
            tmp2 = 0
            for j in range(0, old_data.shape[1]):
                for i in range(0, old_data.shape[0]):
                    old_data_tmp[tmp2].append(old_data[i][j])
                old_data_tmp[tmp2] = np.array(old_data_tmp[tmp2]).reshape(old_data.shape[0],)
                tmp2 += 1
            old_data = np.array(old_data_tmp)
            new_data[0] = new_data[0][0]
            data = []

            for tmp3 in range(0, old_data.shape[0]):
                data_tmp = np.concatenate((old_data[tmp3], new_data[tmp3]), axis=0)
                data.append(data_tmp)
            data = np.array(data)

        elif os.path.isfile(old_monitoring_file) and not os.path.isfile(new_monitoring_file):
            component_count = len(get_configs("velocities",
                                              source + "_" +
                                              get_args("line")).replace(" ", "").split(","))
            components = [i for i in range(1, component_count + 1)]
            old_data = np.loadtxt(old_monitoring_file, dtype=str).reshape(
                (file_len(old_monitoring_file), component_count + 1))
            old_x = correct_numpy_read_data(old_data[:, [0]])
            old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
            old_data[:, [0]] = old_x[0]
            data = old_data

        elif not os.path.isfile(old_monitoring_file) and os.path.isfile(new_monitoring_file):
            component_count = len(get_configs("velocities",
                                              source + "_" +
                                              get_args("line")).replace(" ", "").split(","))
            components = [i for i in range(1, component_count + 1)]
            new_data = np.load(new_monitoring_file, allow_pickle=True)
            data = new_data
        else:
            components = []
            data = None

        for component in components:
            if data is not None:
                index = components.index(component)
                y = data[index + 1]
                y = [np.float128(yi) for yi in y]
                N = len(y)

                variability_index = ((np.max(y) - np.std(y)) - (np.min(y) + np.std(y))) \
                                                  / ((np.max(y) - np.std(y)) + (np.min(y) + np.std(y)))
                function_index = np.sqrt((N / reduce(lambda x, y: x + y, [(1.5 + 0.05 * i) ** 2 for i in y])) *
                                                    ((reduce(lambda x, y: x + y,
                                                              [i ** 2 * (1.5 + 0.05 * i) ** 2 for i in y]) - np.mean(y) *
                                                      reduce(lambda x, y: x + y, [i * (1.5 + 0.05 * i) ** 2 for i in y]))
                                                     / (N - 1)) - 1) / np.mean(y)
                variability_indexes.append(variability_index)
                function_indexes.append(function_index)
                mean_of_y.append(np.mean(y))

    plt.scatter(function_indexes, variability_indexes, c=mean_of_y, cmap="jet")
    plt.ylabel("Variability index")
    plt.xlabel("Function indexes")
    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    main()
