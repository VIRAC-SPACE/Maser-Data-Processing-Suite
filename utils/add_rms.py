import json
import os
import sys

import h5py
import numpy as np

PACKAGE_PARENT = '..'

SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser


def get_configs(section, key):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :return: configuration file section key value
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def split_data_to_signal_and_noise(velocity, amplitude, cuts):
    """

    :param velocity: velocity
    :param amplitude: amplitude
    :param cuts: signal region
    :return: list of signal and list of noise
    """

    cuts_index = list()
    cuts_index.append(0)
    cuts_index2 = list()

    for cut in cuts:
        cuts_index.append((np.abs(velocity - float(cut[0]))).argmin())
        cuts_index.append((np.abs(velocity - float(cut[1]))).argmin())
        cuts_index2.append((np.abs(velocity - float(cut[0]))).argmin())
        cuts_index2.append((np.abs(velocity - float(cut[1]))).argmin())

    cuts_index = sorted(cuts_index)
    cuts_index.append(-1)
    cuts_index2 = sorted(cuts_index2)

    non_signal_amplitude = list()
    signal_amplitude = list()

    i = 0
    j = 1
    while i != len(cuts_index):
        non_signal_amplitude.extend(amplitude[cuts_index[i]: cuts_index[j]])
        i = i + 2
        j = j + 2

    m = 0
    n = 1
    while m != len(cuts_index2):
        signal_amplitude.extend(amplitude[cuts_index2[m]: cuts_index2[n]])
        m = m + 2
        n = n + 2

    return non_signal_amplitude, signal_amplitude


def rms(noise):
    std = np.std(noise)
    number_of_points = len(noise)
    rms = np.sqrt(sum([(n - std) ** 2 for n in noise]) / number_of_points)
    return rms


def main():
    result_files_dir = get_configs("paths", "resultFilePath")
    result_files = os.listdir(result_files_dir)

    for result_file in result_files:
        with open(result_files_dir + result_file) as results:
            if os.stat(result_files_dir + result_file).st_size != 0:
                print(result_file)

                source = result_file.split("_")[0]
                line = result_file.split("_")[1].split(".")[0]
                results_data = json.load(results)
                output_file_path = get_configs("paths", "outputFilePath") + line + "/" + source + "/"
                for experiment in results_data:
                    mjd = results_data[experiment]["modifiedJulianDays"]

                    experiment_data = experiment.split("_")

                    if len(experiment_data) == 7:
                        station = experiment.split("_")[5]
                        index = experiment.split("_")[6]
                    elif len(experiment_data) == 4:
                        station = experiment.split("_")[2]
                        index = experiment.split("_")[3]

                    output_file_name = source + "_" + str(mjd) + "_" + station + "_" + index + ".h5"
                    print(output_file_name)

                    try:
                        output_file = h5py.File(output_file_path + output_file_name, "r")
                        if "amplitude_corrected" in output_file:
                            output_data = output_file["amplitude_corrected"][()]

                            xdata = output_data[:, 0]
                            z1_not_smooht_data = output_data[:, 1]
                            z2_not_smooht_data = output_data[:, 2]
                            avg_y_not_smoohtData = output_data[:, 3]

                            cuts = get_configs('cuts', source + "_" + str(line)).split(";")
                            cuts = [c.split(",") for c in cuts]

                            non_signal_amplitude_left, _ = split_data_to_signal_and_noise(xdata, z1_not_smooht_data,
                                                                                          cuts)
                            non_signal_amplitude_right, _ = split_data_to_signal_and_noise(xdata, z2_not_smooht_data,
                                                                                           cuts)
                            non_signal_amplitude_avg, _ = split_data_to_signal_and_noise(xdata, avg_y_not_smoohtData,
                                                                                         cuts)

                            rms_left = rms(non_signal_amplitude_left)
                            rms_right = rms(non_signal_amplitude_right)
                            rms_avg = rms(non_signal_amplitude_avg)

                            results_data[experiment]["rms_left"] = rms_left
                            results_data[experiment]["rms_right"] = rms_right
                            results_data[experiment]["rms_avg"] = rms_avg

                        else:
                            print("amplitude_corrected not in output file " + output_file_name)
                            results_data[experiment]["rms_left"] = 1.5
                            results_data[experiment]["rms_right"] = 1.5
                            results_data[experiment]["rms_avg"] = 1.5

                    except FileNotFoundError:
                        results_data[experiment]["rms_left"] = 1.5
                        results_data[experiment]["rms_right"] = 1.5
                        results_data[experiment]["rms_avg"] = 1.5
                        print("This file do not exist " + output_file_name)

                with open(result_files_dir + result_file, "w") as output:
                    output.write(json.dumps(results_data, indent=2))


if __name__ == "__main__":
    main()
    sys.exit(0)
