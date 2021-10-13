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


def main(infile):
    plt.style.use('../config/plot.style')
    data = ascii.read(infile)
    sources = data["Source"]
    distance = data["Distance"]
    x = data["X"]
    y = data["Y"]

    Maser = namedtuple('Maser', 'name  x y flux distance '
                                'fluctuation_indexes variability_indexes mean_of_y absolute_mean_of_y, '
                                'largest_y_mean_index outlier')

    masers = [Maser(sources[i], x[i], y[i], [], distance[i], [], [], [], [],
                    [-1], outlier=False) for i in range(0, len(sources))]
    old_monitoring_file_path = get_configs("paths", "oldMonitoringFilePath")
    new_monitoring_file_path = get_configs("paths", "monitoringFilePath")

    for maser in masers:
        print("Executing for maser ", maser.name)

        old = False
        new = False
        both = False

        new_monitoring_file = new_monitoring_file_path + "/" + maser.name + "_" + "6668" + ".npy"
        old_monitoring_file = old_monitoring_file_path + "/" + maser.name + ".dat"

        new_data = None
        old_data = None
        monitoring_data = None
        if os.path.isfile(old_monitoring_file) and os.path.isfile(new_monitoring_file):
            source_velocities = get_configs('velocities', maser.name + "_" + "6668").split(",")
            source_velocities = [si.strip() for si in source_velocities]
            component_count = len(
                get_configs("velocities", maser.name + "_" + "6668").replace(" ", "").split(","))
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
            source_velocities = get_configs('velocities', maser.name + "_" + "6668").split(",")
            source_velocities = [si.strip() for si in source_velocities]
            component_count = len(
                get_configs("velocities", maser.name + "_" + "6668").replace(" ", "").split(","))
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
            source_velocities = get_configs('velocities', maser.name + "_" + "6668").split(",")
            source_velocities = [si.strip() for si in source_velocities]
            component_count = len(
                get_configs("velocities", maser.name + "_" + "6668").replace(" ", "").split(","))
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
                largest_y_mean_index = -1
                largest_y_mean = 0
                for component in components:
                    index = components.index(component)
                    if old:
                        y_data = monitoring_data[:, index + 1]
                    elif both:
                        y_data = monitoring_data[index + 1, :]
                    else:
                        y_data = monitoring_data[index + 1]
                    y_data = np.array([np.float128(yi) for yi in y_data]).clip(min=0)

                    if np.mean(y_data) > largest_y_mean:
                        largest_y_mean = np.mean(y_data)
                        largest_y_mean_index = index
                        maser.largest_y_mean_index[0] = index

                    N = len(y_data)
                    variability_index = ((np.max(y_data) - np.std(y_data)) -
                                         (np.min(y_data) + np.std(y_data))) / \
                                        ((np.max(y_data) - np.std(y_data)) +
                                         (np.min(y_data) + np.std(y_data)))

                    fluctuation_index = np.sqrt(np.abs((N / reduce(lambda x__, y__: x__ + y__,
                                                            [(1.5 + 0.05 * i) ** 2 for i in y_data])) *
                                                ((reduce(lambda x__, y__: x__ + y__,
                                                         [i ** 2 * (1.5 + 0.05 * i) ** 2 for i in y_data]) -
                                                  np.mean(y_data) * reduce(lambda x__, y__: x__ + y__,
                                                                           [i * (1.5 + 0.05 * i) ** 2 for i
                                                                            in y_data])) / (N - 1)) - 1)) / \
                                        np.mean(y_data)

                    if not np.isnan(fluctuation_index):
                        maser.variability_indexes.append(np.float64(variability_index))
                        maser.fluctuation_indexes.append(np.float64(fluctuation_index))
                        maser.mean_of_y.append(np.float64(np.mean(y_data)))
                        maser.flux.extend(y_data)
                    del y_data, variability_index, fluctuation_index

                for component in components:
                    index2 = components.index(component)
                    if index2 == largest_y_mean_index:
                        if old:
                            y_data2 = monitoring_data[:, index2 + 1]
                        elif both:
                            y_data2 = monitoring_data[index2 + 1, :]
                        else:
                            y_data2 = monitoring_data[index2 + 1]
                        y_data2 = np.array([np.float128(yi) for yi in y_data2]).clip(min=0)

                        if maser.distance != "*":
                            maser.absolute_mean_of_y.append(np.float64(np.mean(y_data2)) *
                                                            (np.float64(maser.distance) / 2) ** 2)
                        del y_data2

                del monitoring_data
            del new_monitoring_file, old_monitoring_file

    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig3, ax3 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig4, ax4 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig5, ax5 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig6, ax6 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig7, ax7 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig8, ax8 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig9, ax9 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    fig10, ax10 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)

    x = []
    y = []

    size1 = []
    size2 = []
    size3 = []
    size4 = []
    size5 = []
    size6 = []
    size7 = []

    color5 = []
    color6 = []
    color7 = []

    vis = []
    vis2 = []
    fis = []
    fis2 = []

    distances = []
    absolute_mean_of_ys = []
    mean_of_ys = []
    mean_of_ys2 = []
    for maser in masers:
        if maser.x != "*" and maser.y != "*" and len(maser.mean_of_y) != 0:
            if len(maser.mean_of_y) == 0:
                print("yes", maser.short_name, maser.long_name)
            x.append(float(maser.x))
            y.append(float(maser.y))
            size1.append(np.mean(maser.mean_of_y))
            size2.append(np.mean(maser.absolute_mean_of_y))
            size3.append(100 * np.array(np.mean(maser.fluctuation_indexes)))
            size4.append(100 * np.array(np.mean(maser.variability_indexes)))

        means_y = maser.mean_of_y
        absolute_mean_of_y = maser.absolute_mean_of_y

        for vi in maser.variability_indexes:
            if vi < 0.5:
                color5.append("blue")
            else:
                color5.append("red")

        if len(means_y) > 0:
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
                size5.append(ss)

            vis.extend(maser.variability_indexes)
            fis.extend(maser.fluctuation_indexes)
            mean_of_ys2.extend(maser.mean_of_y)

        if len(absolute_mean_of_y) > 0:
            if 0.5 < absolute_mean_of_y[0] <= 20:
                sss = 100
            elif 20 < absolute_mean_of_y[0] <= 200:
                sss = 200
            elif 200 < absolute_mean_of_y[0] <= 800:
                sss = 300
            elif 800 < absolute_mean_of_y[0] <= 2000:
                sss = 400
            elif 2000 < absolute_mean_of_y[0] <= 3005:
                sss = 500
            else:
                sss = 600
            size6.extend([sss] * len(maser.variability_indexes))
            for vi in maser.variability_indexes:
                if vi < 0.5:
                    color6.append("blue")
                else:
                    color6.append("red")
            vis2.extend(maser.variability_indexes)
            fis2.extend(maser.fluctuation_indexes)

            if maser.distance != "*" and maser.largest_y_mean_index[0] != -1:
                distances.append(float(maser.distance))
                largest_y_mean_index = maser.largest_y_mean_index[0]
                vi_tmp = maser.variability_indexes[largest_y_mean_index]
                if vi_tmp < 0.5:
                    color7.append("blue")
                else:
                    color7.append("red")
                size7.append(sss)
                absolute_mean_of_ys.append(maser.absolute_mean_of_y[0])
                mean_of_ys.append(np.mean(maser.mean_of_y))

    titles_x_y = ["Flux", "CFlux", "Fluctuation indexes", "Variability indexes"]
    ax_x_y = [ax1, ax2, ax3, ax4]
    sizes_x_y = [size1, size2, size3, size4]
    for ax in ax_x_y:
        ax_index = ax_x_y.index(ax)
        ax.scatter(x, y, s=sizes_x_y[ax_index], alpha=0.3)
        ax.set_title(titles_x_y[ax_index])
        ax.set_xlabel("X [Kpc]")
        ax.set_ylabel("Y [Kpc]")

    titles_vi_fi = ["Variability indexes vs Fluctuation indexes Flux",
                    "Variability indexes vs Fluctuation indexes Absolute Flux"]
    ax_vi_fi = [ax5, ax7]
    sizes_vi_fi = [size5, size6]
    colors_vi_fi = [color5, color6]
    viss = [vis, vis2]
    fiss = [fis, fis2]

    for ax in ax_vi_fi:
        ax_index = ax_vi_fi.index(ax)
        scatter = ax.scatter(viss[ax_index], fiss[ax_index], c=colors_vi_fi[ax_index], s=sizes_vi_fi[ax_index],
                             alpha=0.3)
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
        if len(labels) == 5:
            ranges = ["0.5 < Jy <= 20", "20 < Jy <= 200", "200 < Jy <= 800", "800 < Jy <= 2000", "2000 < Jy <= 3005"]
        elif len(labels) == 6:
            ranges = ["0.5 < Jy <= 20", "20 < Jy <= 200", "200 < Jy <= 800", "800 < Jy <= 2000", "2000 < Jy <= 3005",
                      "Jy > 3005"]
        labels = [labels[l] + "  " + ranges[l] for l in range(0, len(ranges))]

        ax.legend(handles, labels, markerscale=0.7)
        ax.set_xlabel("Variability indexes")
        ax.set_ylabel("Fluctuation indexes")
        ax.set_title(titles_vi_fi[ax_index])

    ax_distances_flux_density = [ax6, ax8]
    titles_distances_flux_density = ["Flux vs Distance Absolute Flux", "Flux vs Distance Flux"]
    ax_distances_flux_density_y = [absolute_mean_of_ys, mean_of_ys]
    for ax in ax_distances_flux_density:
        ax_index = ax_distances_flux_density.index(ax)
        ax.set_xlabel("Distance [Kpc]")
        ax.set_ylabel("Flux density (Jy)")
        ax.set_title(titles_distances_flux_density[ax_index])
        ax.scatter(distances, ax_distances_flux_density_y[ax_index], s=size7, c=color7, alpha=0.3)

    ax9.scatter(mean_of_ys2, vis, alpha=0.3, c=color5)
    ax10.scatter(mean_of_ys2, fis, alpha=0.3, c=color5)
    ax9.set_xlabel("Flux density (Jy)")
    ax10.set_xlabel("Flux density (Jy)")
    ax9.set_ylabel("Variability indexes")
    ax10.set_ylabel("Fluctuation indexes")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize parallax data')
    parser.add_argument('inFILE', type=str, help='Specify the input csv file.')
    args = parser.parse_args()
    main(args.inFILE)
    sys.exit(0)
