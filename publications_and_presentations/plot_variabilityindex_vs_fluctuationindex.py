import json
import sys
import os
import argparse
from functools import reduce
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.stats import linregress
import pprint

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser
from utils.help import file_len, convert_datetime_object_to_mjd, find_nearest_index


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
    rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=20)
    new_monitoring_file_path = get_configs("paths", "monitoringFilePath")

    all_sources = list(get_configs_items("sources").keys())

    fluctuation_indexes = []
    variability_indexes = []
    mean_of_y = []
    outliers = []

    for source in all_sources:
        old_monitoring_file = get_configs("paths", "oldMonitoringFilePath") + source + ".dat"
        new_monitoring_file = new_monitoring_file_path + "/" + source + "_" + get_args("line") + ".npy"
        result_file_name = source + "_" + get_args("line") + ".json"
        result_file_path = get_configs("paths", "resultFilePath")
        old_x = []
        old_data = []

        if os.path.isfile(result_file_path + result_file_name):
            if os.path.isfile(new_monitoring_file) and os.path.isfile(old_monitoring_file):
                old_data = np.loadtxt(old_monitoring_file, unpack=True, dtype=str)
                old_x = old_data[0, :]
                old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
                old_data[0, :] = old_x
                new_data = np.load(new_monitoring_file, allow_pickle=True)
                new_x = new_data[0][0]
                x = []
                x.extend(old_x)
                x.extend(new_x)

            elif os.path.isfile(new_monitoring_file) and not os.path.isfile(old_monitoring_file):
                new_data = np.load(new_monitoring_file, allow_pickle=True)
                new_x = new_data[0][0]
                x = new_x

            elif not os.path.isfile(new_monitoring_file) and os.path.isfile(old_monitoring_file):
                old_data = np.loadtxt(old_monitoring_file, unpack=True, dtype=str)
                old_x = old_data[0, :]
                old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
                old_data[0, :] = old_x
                x = old_x

            else:
                continue

            rms = []
            with open(result_file_path + result_file_name) as result_data:
                result = json.load(result_data)
                for experiment in result:
                    rms.append(result[experiment]['rms_avg'])

            if len(old_x) > 0:
                rms_old = [np.mean(rms)] * len(old_x)
                rms_old.extend(rms)
                rms = rms_old

            source_velocities = get_configs('velocities', source + "_" + get_args("line")).split(",")
            source_velocities = [si.strip() for si in source_velocities]
            component_count = len(get_configs("velocities", source + "_" +
                                              get_args("line")).replace(" ", "").split(","))

            components = [i for i in range(1, component_count + 1)]

            if len(old_x) > 0:
                components = components[0:len(old_data) - 1]

            for component in components:
                index = components.index(component)
                if len(old_x) > 0:
                    y_old = old_data[index + 1]
                    y_new = new_data[index + 1]
                    y = []
                    y.extend(y_old)
                    y.extend(y_new)
                else:
                    y_new = new_data[index + 1]
                    y = y_new

                y = np.array([np.float128(yi) for yi in y]).clip(min=0)
                N = len(y)

                error = []
                factor = []
                correction_file = get_configs("parameters", "amplitude_correction_file") + ".npy"
                correction_data = np.load(correction_file)
                correction_mjd = correction_data[:, 0]
                correction_factor = correction_data[:, 1]

                for m in range(0, len(x)):
                    correction_index = find_nearest_index(correction_mjd, x[m])
                    factor.append(correction_factor[correction_index])

                for i in range(0, len(y)):
                    if 1 - factor[i] < 0.05:
                        error.append(2 * rms[i] + y[i] * 0.05)
                    else:
                        error.append(2 * rms[i] + y[i] * (1 - factor[i]))

                variability_index = ((np.max(y) - error[list(y).index(np.max(y))]) - (
                            np.min(y) + error[list(y).index(np.min(y))])) \
                                               / ((np.max(y) - error[list(y).index(np.max(y))]) + (
                            np.min(y) + error[list(y).index(np.min(y))]))

                fluctuation_index = np.sqrt(
                    np.abs((N / reduce(lambda x_, y_: x_ + y_, [error[list(y).index(i)] ** 2 for i in y])) *
                           ((reduce(lambda x_, y_: x_ + y_,
                                    [i ** 2 * (error[list(y).index(i)]) ** 2 for i in y]) -
                             np.mean(y) * reduce(lambda x_, y_: x_ + y_,
                                                 [i * error[list(y).index(i)] ** 2 for i in y]))
                            / (N - 1)) - 1)) / np.mean(y)

                xhi_sqer_red = reduce(lambda x_, y_: x_ + y_,
                                                 [((i - np.mean(y)) / error[list(y).index(i)]) ** 2 for i in y]) / (
                                                      N - 1)

                if not np.isnan(fluctuation_index):
                    variability_indexes.append(np.float64(variability_index))
                    fluctuation_indexes.append(np.float64(fluctuation_index))
                    mean_of_y.append(np.float64(np.mean(y)))
                    if fluctuation_index > 1:
                        source_name = get_configs("Full_source_name", source)
                        outliers.append([variability_index, fluctuation_index, source_name + " " +
                                         source_velocities[index]])

    color = []
    for vi in variability_indexes:
        if vi < 0.5:
            color.append("blue")
        else:
            color.append("red")

    size = []
    for my in mean_of_y:
        if 0.5 < my <= 20:
            size.append(100)
        elif 20 < my <= 200:
            size.append(200)
        elif 200 < my <= 800:
            size.append(300)
        elif 800 < my <= 2000:
            size.append(400)
        elif 2000 < my <= 3005:
            size.append(500)
        else:
            size.append(600)

    scatter = plt.scatter(variability_indexes, fluctuation_indexes, s=size, c=color, alpha=0.3)
    coef = linregress(variability_indexes, fluctuation_indexes)
    plt.plot(variability_indexes, np.array(variability_indexes) * coef.slope + coef.intercept, '--k')
    pprint.pprint(coef, indent=4)

    handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
    ranges = ["0.5 < Jy <= 20", "20 < Jy <= 200", "200 < Jy <= 800", "800 < Jy <= 3000"]
    labels = [labels[l] + "  " + ranges[l] for l in range(0, len(ranges))]

    plt.legend(handles, labels)
    for outlier in outliers:
        plt.text(outlier[0] - 0.045, outlier[1] + 0.045, outlier[2], fontsize=11)

    plt.xlabel("Variability index")
    plt.ylabel("Fluctuation index")
    plt.show()

    sys.exit(0)


if __name__ == "__main__":
    main()
