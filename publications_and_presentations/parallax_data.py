import sys
import os
from collections import namedtuple
from functools import reduce
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


def get_maser_short_name(source):
    all_sources_short_name = list(get_configs_items("Full_source_name").keys())
    all_sources_full_name = list(get_configs_items("Full_source_name").values())
    try:
        source_short_name = all_sources_short_name[all_sources_full_name.index(source)]
    except ValueError:
        source_short_name = "*"
        print(source, " is not in config file")

    return source_short_name


def main(infile):
    data = ascii.read(infile)
    sources = data["Source"]
    distance = data["Distance"]
    x = data["X"]
    y = data["Y"]

    Maser = namedtuple('Maser', 'long_name short_name x y flux distance '
                                'fluctuation_indexes variability_indexes mean_of_y absolute_mean_of_y, outlier')

    masers = [Maser(sources[i], get_maser_short_name(sources[i]), x[i], y[i], [], distance[i], [], [], [], [],
                    outlier=False) for i in range(0, len(sources))]

    old_monitoring_file_path = get_configs("paths", "oldMonitoringFilePath")
    new_monitoring_file_path = get_configs("paths", "monitoringFilePath")

    for maser in masers:
        print("Executing for maser with short name", maser.short_name, "long name", maser.long_name)

        if maser.short_name != "*":
            old = False
            new = False
            both = False

            new_monitoring_file = new_monitoring_file_path + "/" + maser.short_name + "_" + "6668" + ".npy"
            old_monitoring_file = old_monitoring_file_path + "/" + maser.short_name + ".dat"

            new_data = None
            old_data = None
            monitoring_data = None
            if os.path.isfile(old_monitoring_file) and os.path.isfile(new_monitoring_file):
                source_velocities = get_configs('velocities', maser.short_name + "_" + "6668").split(",")
                source_velocities = [si.strip() for si in source_velocities]
                component_count = len(
                    get_configs("velocities", maser.short_name + "_" + "6668").replace(" ", "").split(","))
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
                source_velocities = get_configs('velocities', maser.short_name + "_" + "6668").split(",")
                source_velocities = [si.strip() for si in source_velocities]
                component_count = len(
                    get_configs("velocities", maser.short_name + "_" + "6668").replace(" ", "").split(","))
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
                source_velocities = get_configs('velocities', maser.short_name + "_" + "6668").split(",")
                source_velocities = [si.strip() for si in source_velocities]
                component_count = len(
                    get_configs("velocities", maser.short_name + "_" + "6668").replace(" ", "").split(","))
                components = [i for i in range(1, component_count + 1)]
                new_data = np.load(new_monitoring_file, allow_pickle=True)
                monitoring_data = new_data
                new = True

            else:
                components = []
                source_velocities = []
                data = None

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
                        variability_index = ((np.max(y_data) - np.std(y_data)) -
                                             (np.min(y_data) + np.std(y_data))) / \
                                            ((np.max(y_data) - np.std(y_data)) +
                                             (np.min(y_data) + np.std(y_data)))

                        fluctuation_index = np.sqrt((N / reduce(lambda x__, y__: x__ + y__,
                                                                [(1.5 + 0.05 * i) ** 2 for i in y_data])) *
                                                    ((reduce(lambda x__, y__: x__ + y__,
                                                             [i ** 2 * (1.5 + 0.05 * i) ** 2 for i in y_data]) -
                                                      np.mean(y_data) * reduce(lambda x__, y__: x__ + y__,
                                                                               [i * (1.5 + 0.05 * i) ** 2 for i
                                                                                in y_data])) / (N - 1)) - 1) / \
                                            np.mean(y_data)
                        if not np.isnan(fluctuation_index):
                            maser.variability_indexes.append(np.float64(variability_index))
                            maser.fluctuation_indexes.append(np.float64(fluctuation_index))
                            maser.mean_of_y.append(np.float64(np.mean(y_data)))
                            maser.flux.extend(y_data)
                            if maser.distance != "*":
                                maser.absolute_mean_of_y.append(np.float64(np.mean(y_data)) *
                                                                (np.float64(maser.distance) / 2) ** 2)

                            #if fluctuation_index > 1:
                                #maser.outlier = True

                        del y_data, variability_index, fluctuation_index
                    del monitoring_data
                del new_monitoring_file, old_monitoring_file

    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig3, ax3 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig4, ax4 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig5, ax5 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig6, ax6 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig7, ax7 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)

    axs = [ax1, ax2, ax3, ax4]
    for ax in axs:
        ax.set_xlabel("X [Kpc]")
        ax.set_ylabel("Y [Kpc]")

    ax5.set_xlabel("Variability indexes")
    ax5.set_ylabel("Fluctuation indexes")
    ax5.set_title("Variability indexes vs Fluctuation indexes")

    ax6.set_xlabel("Distance [Kpc]")
    ax6.set_ylabel("Flux density (Jy)")
    ax6.set_title("Flux vs Distance")

    ax7.xlabel("Variability indexes")
    ax7.ylabel("Fluctuation indexes")
    ax7.title("Variability indexes vs Fluctuation indexes")

    for maser in masers:
        if maser.x != "*" and maser.y != "*":
            ax1.scatter(float(maser.x), float(maser.y), s=np.mean(maser.mean_of_y))
            ax1.set_title("Flux")

            ax2.scatter(float(maser.x), float(maser.y), s=np.mean(maser.absolute_mean_of_y))
            ax2.set_title("CFlux")

            ax3.scatter(float(maser.x), float(maser.y), s=100 * np.array(np.mean(maser.fluctuation_indexes)))
            ax3.set_title("Fluctuation indexes")

            ax4.scatter(float(maser.x), float(maser.y), s=100 * np.array(np.mean(maser.variability_indexes)))
            ax4.set_title("Variability indexes")

        collor1 = []
        size1 = []
        size2 = []
        variability_indexes = maser.variability_indexes
        means_y = maser.mean_of_y
        absolute_mean_of_y = maser.absolute_mean_of_y
        for vi in variability_indexes:
            if vi < 0.5:
                collor1.append("blue")
            else:
                collor1.append("red")

        for my in means_y:
            if 0.5 < my <= 20:
                ss = 100
            elif 20 < my <= 200:
                ss = 200
            elif 200 < my <= 800:
               ss = 300
            elif 800 < my <= 2000:
                ss = 400
            elif 2000 < my <= 3005:
                ss = 500
            else:
                ss = 600
            size1.append(ss)

            if maser.distance != "*":
                my_index = means_y.index(my)
                vi_tmp = variability_indexes[my_index]
                if vi_tmp < 0.5:
                    ax6.scatter(float(maser.distance), my, s=ss, c="b")
                else:
                    ax6.scatter(float(maser.distance), my, s=ss, c="r")

        for my2 in absolute_mean_of_y:
            if 0.5 < my2 <= 20:
                sss = 100
            elif 20 < my2 <= 200:
                sss = 200
            elif 200 < my2 <= 800:
               sss = 300
            elif 800 < my2 <= 2000:
                sss = 400
            elif 2000 < my2 <= 3005:
                sss = 500
            else:
                sss = 600
            size2.append(sss)

        ax5.scatter(maser.variability_indexes, maser.fluctuation_indexes, c=collor1, s=size1)
        ax7.scatter(maser.variability_indexes, maser.fluctuation_indexes, c=collor1, s=size2)

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize parallax data')
    parser.add_argument('inFILE', type=str, help='Specify the input csv file.')
    args = parser.parse_args()
    main(args.inFILE)
    sys.exit(0)
