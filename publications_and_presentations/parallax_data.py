import json
import sys
import os
from collections import namedtuple

from functools import reduce
import argparse
from datetime import datetime
from astropy.io import ascii
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import FK5, galactocentric_frame_defaults
import numpy as np
import matplotlib.pyplot as plt
from mw_plot import MWPlot

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser
from utils.help import file_len, correct_numpy_read_data, convert_datetime_object_to_mjd, find_nearest_index


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
    galactocentric_frame_defaults.set('latest')
    plt.style.use('../config/plot.style')
    data = ascii.read(infile)
    sources = data["Source"]
    distance = data["Distance"]
    x = data["X"]
    y = data["Y"]

    Maser = namedtuple('Maser', 'name  x y flux distance '
                                'fluctuation_indexes variability_indexes mean_of_y absolute_mean_of_y, '
                                'largest_y_mean_index  mean_of_y2 xhi_sqer_red outlier')

    masers = [Maser(sources[i], x[i], y[i], [], distance[i], [], [], [], [],
                    [-1], [], [], outlier=False) for i in range(0, len(sources))]
    old_monitoring_file_path = get_configs("paths", "oldMonitoringFilePath")
    new_monitoring_file_path = get_configs("paths", "monitoringFilePath")

    for maser in masers:
        print("Executing for maser ", maser.name)

        old = False
        new = False
        both = False

        new_monitoring_file = new_monitoring_file_path + "/" + maser.name + "_" + "6668" + ".npy"
        old_monitoring_file = old_monitoring_file_path + "/" + maser.name + ".dat"

        result_file_name = maser.name + "_" + "6668" + ".json"
        result_file_path = get_configs("paths", "resultFilePath")

        new_data = None
        old_data = None
        monitoring_data = None
        old_x = []
        if os.path.isfile(result_file_path + result_file_name):
            if os.path.isfile(new_monitoring_file) and os.path.isfile(old_monitoring_file):
                old_data = np.loadtxt(old_monitoring_file, unpack=True, dtype=str)
                old_x = old_data[0, :]
                old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
                old_data[0, :] = old_x
                new_data = np.load(new_monitoring_file, allow_pickle=True)
                new_x = new_data[0][0]
                x_mjd = []
                x_mjd.extend(old_x)
                x_mjd.extend(new_x)

                monitoring_data = []
                data_tmp = np.concatenate((old_data[0, :], new_data[0][0]), axis=0)
                monitoring_data.append(data_tmp)

                for tmp3 in range(1, old_data.shape[0]):
                    data_tmp = np.concatenate((old_data[tmp3, :], new_data[tmp3]), axis=0)
                    monitoring_data.append(data_tmp)
                monitoring_data = np.array(monitoring_data)
                both = True

            elif os.path.isfile(new_monitoring_file) and not os.path.isfile(old_monitoring_file):
                new_data = np.load(new_monitoring_file, allow_pickle=True)
                new_x = new_data[0][0]
                x_mjd = new_x
                monitoring_data = new_data
                new = True

            elif not os.path.isfile(new_monitoring_file) and os.path.isfile(old_monitoring_file):
                old_data = np.loadtxt(old_monitoring_file, unpack=True, dtype=str)
                old_x = old_data[0, :]
                old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
                old_data[0, :] = old_x
                x_mjd = old_x
                monitoring_data = old_data
                old = True

            else:
                monitoring_data = None
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

            source_velocities = get_configs('velocities', maser.name + "_" + "6668").split(",")
            source_velocities = [si.strip() for si in source_velocities]
            component_count = len(get_configs("velocities", maser.name + "_" + "6668").replace(" ", "").split(","))
            components = [i for i in range(1, component_count + 1)]

            if len(old_x) > 0:
                components = components[0:len(old_data) - 1]

            if len(components) > 0:
                if monitoring_data is not None:
                    for component in components:

                        index = components.index(component)
                        if len(old_x) > 0:
                            y_old = old_data[index + 1]
                            y_new = new_data[index + 1]
                            y_data = []
                            y_data.extend(y_old)
                            y_data.extend(y_new)
                        else:
                            y_new = new_data[index + 1]
                            y_data = y_new

                        y_data = np.array([np.float128(yi) for yi in y_data]).clip(min=0)

                        largest_y_mean_index = -1
                        largest_y_mean = 0

                        if np.mean(y_data) > largest_y_mean:
                            largest_y_mean = np.mean(y_data)
                            largest_y_mean_index = index
                            maser.largest_y_mean_index[0] = index

                        N = len(y_data)
                        rms = []
                        with open(result_file_path + result_file_name) as result_data:
                            result = json.load(result_data)
                            for experiment in result:
                                rms.append(result[experiment]['rms_avg'])

                        if len(old_x) > 0:
                            rms_old = [np.mean(rms)] * len(old_x)
                            rms_old.extend(rms)
                            rms = rms_old

                        error = []
                        factor = []
                        correction_file = get_configs("parameters", "amplitude_correction_file") + ".npy"
                        correction_data = np.load(correction_file)
                        correction_mjd = correction_data[:, 0]
                        correction_factor = correction_data[:, 1]

                        for m in range(0, len(x_mjd)):
                            correction_index = find_nearest_index(correction_mjd, x_mjd[m])
                            factor.append(correction_factor[correction_index])

                        for i in range(0, len(y_data)):
                            if 1 - factor[i] < 0.05:
                                error.append(2 * rms[i] + y_data[i] * 0.05)
                            else:
                                error.append(2 * rms[i] + y_data[i] * (1 - factor[i]))

                        variability_index = ((np.max(y_data) - error[list(y_data).index(np.max(y_data))]) - (
                                np.min(y_data) + error[list(y_data).index(np.min(y_data))])) \
                                            / ((np.max(y_data) - error[list(y_data).index(np.max(y_data))]) + (
                                np.min(y_data) + error[list(y_data).index(np.min(y_data))]))

                        fluctuation_index = np.sqrt(
                            np.abs((N / reduce(lambda x_, y_: x_ + y_,
                                               [error[list(y_data).index(i)] ** 2 for i in y_data])) *
                                   ((reduce(lambda x_, y_: x_ + y_,
                                            [i ** 2 * (error[list(y_data).index(i)]) ** 2 for i in y_data]) -
                                     np.mean(y_data) * reduce(lambda x_, y_: x_ + y_,
                                                         [i * error[list(y_data).index(i)] ** 2 for i in y_data]))
                                    / (N - 1)) - 1)) / np.mean(y_data)

                        xhi_sqer_red = reduce(lambda x_, y_: x_ + y_,
                                              [((i - np.mean(y_data)) /
                                                error[list(y_data).index(i)]) ** 2 for i in y_data]) / (
                                               N - 1)

                        maser.xhi_sqer_red.append(xhi_sqer_red)

                        if not np.isnan(fluctuation_index):
                            maser.variability_indexes.append(np.float64(variability_index))
                            maser.fluctuation_indexes.append(np.float64(fluctuation_index))
                            maser.mean_of_y.append(np.float64(np.mean(y_data)))
                            maser.flux.extend(y_data)
                        del y_data, variability_index, fluctuation_index

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

                            maser.mean_of_y2.append(np.float64(np.mean(y_data2)))
                            del y_data2

    mw1 = MWPlot(rot90=45, r0=8.5, radius=20 * u.kpc, unit=u.kpc, coord='galactocentric', annotation=True,
                 figsize=(16, 16), dpi=150)
    mw2 = MWPlot(rot90=45, r0=8.5, radius=20 * u.kpc, unit=u.kpc, coord='galactocentric', annotation=True,
                 figsize=(16, 16), dpi=150)
    mw3 = MWPlot(rot90=45, r0=8.5, radius=20 * u.kpc, unit=u.kpc, coord='galactocentric', annotation=True,
                 figsize=(16, 16), dpi=150)
    mw4 = MWPlot(rot90=45, r0=8.5, radius=20 * u.kpc, unit=u.kpc, coord='galactocentric', annotation=True,
                 figsize=(16, 16), dpi=150)

    fig5, ax5 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)  # fig1
    fig6, ax6 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)  # fig2
    fig7, ax7 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)  # fig3
    fig8, ax8 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)  # fig4
    fig9, ax9 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)  # fig5
    fig10, ax10 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)  # fig6
    fig11, ax11 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)  # fig7
    fig12, ax12 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)  # fig8
    fig13, ax13 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)  # fig9
    fig14, ax14 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)  # fig10

    x_ = []
    y_ = []

    color_ = []
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
    xhi_sqer = []

    distances = []
    absolute_mean_of_ys = []
    mean_of_ys = []
    mean_of_ys2 = []
    mean_of_ys3 = []

    for maser in masers:
        if maser.distance != "*":
            if np.mean(maser.variability_indexes) < 0.5:
                c = "blue"
            else:
                c = "red"
            if len(maser.mean_of_y) != 0:
                size1.append(np.mean(maser.mean_of_y))
                size2.append(np.mean(maser.absolute_mean_of_y))
                size3.append(100 * np.array(np.mean(maser.fluctuation_indexes)))
                size4.append(100 * np.array(np.mean(maser.variability_indexes)))
                source_cordinations = get_configs("sources", maser.name).split(",")
                source_cordinations = [sc.strip() for sc in source_cordinations]
                RA = source_cordinations[0]
                DEC = source_cordinations[1]

                ra = list()
                dec = list()
                ra.append(RA[0:2])
                ra.append(RA[2:4])
                ra.append(RA[4:len(RA)])

                if DEC[0] == "-":
                    dec.append(DEC[0:3])
                    dec.append(DEC[3:5])
                    dec.append(DEC[5:len(DEC)])
                else:
                    dec.append(DEC[0:2])
                    dec.append(DEC[2:4])
                    dec.append(DEC[4:len(DEC)])

                ra_str = ra[0] + "h" + ra[1] + "m" + ra[2] + "s"
                if int(dec[0]) > 0:
                    dec_str = "+" + dec[0] + "d" + dec[1] + "m" + dec[2] + "s"
                else:
                    dec_str = dec[0] + "d" + dec[1] + "m" + dec[2] + "s"

                coord_ = coord.SkyCoord(ra=ra_str, dec=dec_str, frame=FK5, equinox='J2000.0')
                coord_ = coord_.galactic
                lon = coord_.l.deg
                x_.append(np.sin(np.radians(lon)) * float(maser.distance))
                y_.append(8.5 - np.cos(np.radians(lon)) * float(maser.distance))

                vi = np.mean(maser.variability_indexes)
                if vi < 0.5:
                    color_.append("blue")
                else:
                    color_.append("red")

            ax11.scatter(float(maser.distance), np.mean(maser.fluctuation_indexes), c=c, alpha=0.3)
            ax12.scatter(float(maser.distance), np.mean(maser.variability_indexes), c=c, alpha=0.3)

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
            xhi_sqer.extend(maser.xhi_sqer_red)
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
                mean_of_ys3.append(maser.mean_of_y2[0])

    titles_x_y = ["Flux", "CFlux", "Fluctuation indexes", "Variability indexes"]
    ax_x_y = [mw1, mw2, mw3, mw4]
    sizes_x_y = [size1, size2, size3, size4]

    sun = (0, 8.5)
    gc = (0, 0)
    color = ["black", "goldenrod"] + color_
    x__ = [0.0, 8.5] + y_
    y__ = [0.0, 0.0] + x_
    for ax in ax_x_y:
        ax_index = ax_x_y.index(ax)
        size = [20, 20] + sizes_x_y[ax_index]
        ax.scatter(x=x__ * u.kpc, y=y__ * u.kpc, marker="o", s=size, c=color, alpha=0.3)
        ax.title = titles_x_y[ax_index]
        #ax.set_xlabel("X [Kpc]")
        #ax.set_ylabel("Y [Kpc]")

    titles_vi_fi = ["Variability indexes vs Fluctuation indexes Flux",
                    "Variability indexes vs Fluctuation indexes Absolute Flux"]
    ax_vi_fi = [ax5, ax7]
    sizes_vi_fi = [size5, size6]
    colors_vi_fi = [color5, color6]
    viss = [vis, vis2]
    fiss = [fis, fis2]

    for ax in ax_vi_fi:
        ax_index = ax_vi_fi.index(ax)

        scatter = ax.scatter(viss[ax_index], fiss[ax_index], c=colors_vi_fi[ax_index],
                             s=sizes_vi_fi[ax_index], alpha=0.3)
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
    ax_distances_flux_density_y = [absolute_mean_of_ys, mean_of_ys3]
    for ax in ax_distances_flux_density:
        ax_index = ax_distances_flux_density.index(ax)
        ax.set_xlabel("Distance [Kpc]")
        ax.set_ylabel("Flux density (Jy)")
        ax.set_title(titles_distances_flux_density[ax_index])
        ax.scatter(distances, ax_distances_flux_density_y[ax_index], s=size7, c=color7, alpha=0.3)

    ax9.scatter(mean_of_ys2, vis, alpha=0.3, c=color5)
    ax10.scatter(mean_of_ys2, fis, alpha=0.3, c=color5)

    ax13.scatter(xhi_sqer, vis, alpha=0.3, c=color5)
    ax14.scatter(xhi_sqer, fis, alpha=0.3, c=color5)

    ax9.set_xlabel("Flux density (Jy)")
    ax10.set_xlabel("Flux density (Jy)")
    ax9.set_ylabel("Variability indexes")
    ax10.set_ylabel("Fluctuation indexes")
    ax11.set_xlabel("Distance [Kpc]")
    ax11.set_ylabel("Variability indexes")
    ax12.set_xlabel("Distance [Kpc]")
    ax12.set_ylabel("Fluctuation indexes")
    ax13.set_xlabel("xhi sqer red")
    ax13.set_ylabel("Variability indexes")
    ax14.set_xlabel("xhi sqer red")
    ax14.set_ylabel("Fluctuation indexes")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize parallax data')
    parser.add_argument('inFILE', type=str, help='Specify the input csv file.')
    args = parser.parse_args()
    main(args.inFILE)
    sys.exit(0)
