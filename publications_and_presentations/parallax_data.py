from functools import reduce

import sys
import os
import argparse
from datetime import datetime

from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser
from utils.help import file_len, correct_numpy_read_data, convert_datetime_object_to_mjd


def get_configs(section, key):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :return: configuration file section key value
    """
    config_file_path = "../config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def get_configs_items(section):
    """

    :return: None
    """
    config_file_path = "../config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_items(section)


def main(infile):
    data = ascii.read(infile)
    sources = data["Source"]
    flux = data["Flux"]
    cflux = data["Cflux"]
    distance = data["Distance"]
    x_ = np.array([float(data["X"][i]) for i in range(0, len(data["X"])) if data["X"][i] != "*"])
    y_ = np.array([float(data["Y"][i]) for i in range(0, len(data["Y"])) if data["Y"][i] != "*"])
    flux = np.array([float(flux[i]) for i in range(0, len(flux)) if data["X"][i] != "*"])
    cflux = np.array([float(cflux[i]) for i in range(0, len(cflux)) if data["Y"][i] != "*"])
    sources = [sources[i] for i in range(0, len(sources)) if data["Y"][i] != "*"]
    distance = [float(distance[i]) for i in range(0, len(distance)) if data["Y"][i] != "*"]

    old_monitoring_file_path = get_configs("paths", "oldMonitoringFilePath")
    new_monitoring_file_path = get_configs("paths", "monitoringFilePath")

    all_sources_short_name = list(get_configs_items("Full_source_name").keys())
    all_sources_full_name = list(get_configs_items("Full_source_name").values())

    fluctuation_indexes = dict()
    variability_indexes = dict()
    mean_of_y = dict()
    absolute_mean_of_y = dict()
    distance_ = dict()
    outliers = []

    for source in sources:
        old = False
        new = False
        both = False
        print(source)
        try:
            source_short_name = all_sources_short_name[all_sources_full_name.index(source)]
        except ValueError:
            print(source, " is not in config file")
        else:

            new_monitoring_file = new_monitoring_file_path + "/" + source_short_name + "_" + "6668" + ".npy"
            old_monitoring_file = old_monitoring_file_path + "/" + source_short_name + ".dat"
            fluctuation_indexes[source_short_name] = []
            variability_indexes[source_short_name] = []
            mean_of_y[source_short_name] = []
            absolute_mean_of_y[source_short_name] = []
            distance_[source_short_name] = distance[sources.index(source)]

            new_data = None
            old_data = None
            monitoring_data = None

            if os.path.isfile(old_monitoring_file) and os.path.isfile(new_monitoring_file):
                source_velocities = get_configs('velocities', source_short_name + "_" + "6668").split(",")
                source_velocities = [si.strip() for si in source_velocities]
                component_count = len(get_configs("velocities", source_short_name + "_" + "6668").replace(" ", "").split(","))
                components = [i for i in range(1, component_count + 1)]
                new_data = np.load(new_monitoring_file, allow_pickle=True)
                new_x = new_data[0][0]
                old_data = np.loadtxt(old_monitoring_file, dtype=str).reshape(
                    (file_len(old_monitoring_file), component_count + 1))
                old_x = correct_numpy_read_data(old_data[:, [0]])
                old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
                old_data[:, 0] = old_x
                x_data = old_x + list(new_x)
                old_data_tmp = []
                for tmp in range(0, len(new_data)):
                    old_data_tmp.append([])
                tmp2 = 0
                for j in range(0, old_data.shape[1]):
                    for i in range(0, old_data.shape[0]):
                        old_data_tmp[tmp2].append(old_data[i][j])
                    old_data_tmp[tmp2] = np.array(old_data_tmp[tmp2]).reshape(old_data.shape[0])
                    tmp2 += 1
                old_data = np.array(old_data_tmp)
                new_data[0] = new_data[0][0]
                monitoring_data = []

                for tmp3 in range(0, old_data.shape[0]):
                    data_tmp = np.concatenate((old_data[tmp3], new_data[tmp3]), axis=0)
                    monitoring_data.append(data_tmp)
                monitoring_data = np.array(monitoring_data)
                both = True

            elif os.path.isfile(old_monitoring_file) and not os.path.isfile(new_monitoring_file):
                source_velocities = get_configs('velocities', source_short_name + "_" + "6668").split(",")
                source_velocities = [si.strip() for si in source_velocities]
                component_count = len(get_configs("velocities", source_short_name + "_" + "6668").replace(" ", "").split(","))
                components = [i for i in range(1, component_count + 1)]
                old_data = np.loadtxt(old_monitoring_file, dtype=str).reshape(
                    (file_len(old_monitoring_file), component_count + 1))
                old_x = correct_numpy_read_data(old_data[:, [0]])
                old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
                old_data[:, 0] = old_x
                monitoring_data = old_data
                x = list(old_x)
                old = True

            elif not os.path.isfile(old_monitoring_file) and os.path.isfile(new_monitoring_file):
                source_velocities = get_configs('velocities', source_short_name + "_" + "6668").split(",")
                source_velocities = [si.strip() for si in source_velocities]
                component_count = len(get_configs("velocities", source_short_name + "_" + "6668").replace(" ", "").split(","))
                components = [i for i in range(1, component_count + 1)]
                new_data = np.load(new_monitoring_file, allow_pickle=True)
                monitoring_data = new_data
                new = True

            else:
                components = []
                source_velocities = []
                monitoring_data = None

            if len(components) > 0:
                if monitoring_data is not None:
                    for component in components:
                        index = components.index(component)
                        if old:
                            y_data = monitoring_data[:, index + 1]
                        elif both:
                            y_data = monitoring_data[index + 1, :]
                        else:
                            y_data = monitoring_data[index + 1]
                        y_data = np.array([np.float128(yi) for yi in y_data]).clip(min=0)
                        N = len(y_data)
                        variability_index = ((np.max(y_data) - np.std(y_data)) - (np.min(y_data) + np.std(y_data))) / \
                                            ((np.max(y_data) - np.std(y_data)) + (np.min(y_data) + np.std(y_data)))
                        fluctuation_index = np.sqrt((N / reduce(lambda x__, y__: x__ + y__, [(1.5 + 0.05 * i) ** 2 for i in y_data])) *
                                                    ((reduce(lambda x__, y__: x__ + y__,
                                                             [i ** 2 * (1.5 + 0.05 * i) ** 2 for i in y_data]) - np.mean(y_data) *
                                                      reduce(lambda x__, y__: x__ + y__, [i * (1.5 + 0.05 * i) ** 2 for i in y_data]))
                                                     / (N - 1)) - 1) / np.mean(y_data)
                        if not np.isnan(fluctuation_index):
                            fluctuation_indexes[source_short_name].append(np.float64(variability_index))
                            variability_indexes[source_short_name].append(np.float64(fluctuation_index))
                            mean_of_y[source_short_name].append(np.float64(np.mean(y_data)))
                            distance_tmp = distance[sources.index(source)]
                            absolute_mean_of_y[source_short_name].append(np.float64(np.mean(y_data)) * (distance_tmp/2)**2)
                            if fluctuation_index > 1:
                                source_name = get_configs("Full_source_name", source_short_name)
                                outliers.append(
                                    [variability_index, fluctuation_index, source_name + " " + source_velocities[index]])

                        del y_data, variability_index, fluctuation_index
                    del monitoring_data
            del new_monitoring_file, old_monitoring_file

    fluctuation_indexes_sizes = []
    variability_indexes_sizes = []
    fluctuation_indexes_vs_variability_indexes_sizes = []
    fluctuation_indexes_vs_variability_indexes_sizes_absolute = []
    fluctuation_indexes_vs_variability_color = []
    mean_of_y_ = []

    for source in sources:
        try:
            source_short_name = all_sources_short_name[all_sources_full_name.index(source)]
        except ValueError:
            print(source, " is not in config file")
        else:
            fluctuation_indexes_sizes.append(np.mean(fluctuation_indexes[source_short_name]))
            variability_indexes_sizes.append(np.mean(variability_indexes[source_short_name]))
            mean_of_y_.append(np.mean(mean_of_y[source_short_name]))

            vi = np.mean(variability_indexes[source_short_name])
            if vi < 0.5:
                fluctuation_indexes_vs_variability_color.append("blue")
            else:
                fluctuation_indexes_vs_variability_color.append("red")

            my = np.mean(mean_of_y[source_short_name])
            if 0.5 < my <= 20:
                fluctuation_indexes_vs_variability_indexes_sizes.append(100)
            elif 20 < my <= 200:
                fluctuation_indexes_vs_variability_indexes_sizes.append(200)
            elif 200 < my <= 800:
                fluctuation_indexes_vs_variability_indexes_sizes.append(300)
            elif 800 < my <= 2000:
                fluctuation_indexes_vs_variability_indexes_sizes.append(400)
            elif 2000 < my <= 3005:
                fluctuation_indexes_vs_variability_indexes_sizes.append(500)
            else:
                fluctuation_indexes_vs_variability_indexes_sizes.append(600)

            my = np.mean(absolute_mean_of_y[source_short_name])
            if 0.5 < my <= 20:
                fluctuation_indexes_vs_variability_indexes_sizes_absolute.append(100)
            elif 20 < my <= 200:
                fluctuation_indexes_vs_variability_indexes_sizes_absolute.append(200)
            elif 200 < my <= 800:
                fluctuation_indexes_vs_variability_indexes_sizes_absolute.append(300)
            elif 800 < my <= 2000:
                fluctuation_indexes_vs_variability_indexes_sizes_absolute.append(400)
            elif 2000 < my <= 3005:
                fluctuation_indexes_vs_variability_indexes_sizes_absolute.append(500)
            else:
                fluctuation_indexes_vs_variability_indexes_sizes_absolute.append(600)

    plt.figure(dpi=150)
    plt.scatter(x_, y_, s=0.1 * flux)
    plt.xlim(min(x_) - 0.3, max(x_) + 0.3)
    plt.ylim(min(y_) - 0.3, max(y_) + 0.3)
    plt.xlabel("X [Kpc]")
    plt.ylabel("Y [Kpc]")
    plt.title("Flux")

    plt.figure(dpi=150)
    plt.scatter(x_, y_, s=0.1 * cflux)
    plt.xlim(min(x_) - 0.3, max(x_) + 0.3)
    plt.ylim(min(y_) - 0.3, max(y_) + 0.3)
    plt.xlabel("X [Kpc]")
    plt.ylabel("Y [Kpc]")
    plt.title("CFlux")

    plt.figure(dpi=150)
    plt.scatter(x_, y_, s=100 * np.array(fluctuation_indexes_sizes))
    plt.xlim(min(x_) - 0.3, max(x_) + 0.3)
    plt.ylim(min(y_) - 0.3, max(y_) + 0.3)
    plt.xlabel("X [Kpc]")
    plt.ylabel("Y [Kpc]")
    plt.title("Fluctuation indexes")

    plt.figure(dpi=150)
    plt.scatter(x_, y_, s=100 * np.array(variability_indexes_sizes))
    plt.xlim(min(x_) - 0.3, max(x_) + 0.3)
    plt.ylim(min(y_) - 0.3, max(y_) + 0.3)
    plt.xlabel("X [Kpc]")
    plt.ylabel("Y [Kpc]")
    plt.title("Variability indexes")

    plt.figure(dpi=150)
    plt.scatter(variability_indexes_sizes, fluctuation_indexes_sizes,
                s=fluctuation_indexes_vs_variability_indexes_sizes, c=fluctuation_indexes_vs_variability_color)
    plt.xlabel("Variability indexes")
    plt.ylabel("Fluctuation indexes")
    plt.title("Variability indexes vs Fluctuation indexes")

    plt.figure(dpi=150)
    plt.scatter(distance, mean_of_y_,
                s=fluctuation_indexes_vs_variability_indexes_sizes, c=fluctuation_indexes_vs_variability_color)
    plt.xlabel("Distance [Kpc]")
    plt.ylabel("Flux density (Jy)")
    plt.title("Flux vs Distance")

  
    plt.figure(dpi=150)
    plt.scatter(variability_indexes_sizes, fluctuation_indexes_sizes,
                s=fluctuation_indexes_vs_variability_indexes_sizes_absolute, c=fluctuation_indexes_vs_variability_color)
    plt.xlabel("Variability indexes")
    plt.ylabel("Fluctuation indexes")
    plt.title("Variability indexes vs Fluctuation indexes")

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize parallax data')
    parser.add_argument('inFILE', type=str, help='Specify the input csv file.')
    args = parser.parse_args()
    main(args.inFILE)
    sys.exit(0)
