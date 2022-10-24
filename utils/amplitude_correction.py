import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import h5py


PACKAGE_PARENT = '..'

SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))


from parsers.configparser_ import ConfigParser
from utils.help import find_nearest_index


def get_configs(section, key):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :return: configuration file section key value
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def main():
    correction_file = get_configs("parameters", "amplitude_correction_file") + ".npy"
    correction_data = np.load(correction_file)
    correction_mjd = correction_data[:, 0]
    correction_factor = correction_data[:, 1]

    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 8), dpi=90)
    ax1.scatter(correction_data[:, 0], correction_data[:, 1])
    ax1.set_xlabel("mjd")
    ax1.set_ylabel("amplitude")

    result_files_dir = get_configs("paths", "resultFilePath")
    result_files = os.listdir(result_files_dir)

    for result_file in result_files:
        print(result_file)
        with open(result_files_dir + result_file) as results:
            if os.stat(result_files_dir + result_file).st_size != 0:
                results_data = json.load(results)

                source = result_file.split("_")[0]
                line = result_file.split("_")[1].split(".")[0]

                source_velocities = get_configs('velocities', source + "_" + line).split(",")
                source_velocities = [x.strip() for x in source_velocities]

                output_file_path = get_configs("paths", "outputFilePath") + line + "/" + source + "/"

                for experiment in results_data:
                    mjd = results_data[experiment]["modifiedJulianDays"]
                    correction_index = find_nearest_index(correction_mjd, float(mjd))
                    factor = correction_factor[correction_index]

                    station = experiment.split("_")[5]
                    index = experiment.split("_")[6]
                    output_file_name = source + "_" + str(mjd) + "_" + station + "_" + index + ".h5"
                    print(output_file_name)

                    output_file = h5py.File(output_file_path + output_file_name, "a")
                    if "amplitude_corrected" in output_file and "gain_corrected" not in output_file:
                        gain_corrected_results = output_file["amplitude_corrected"][()]
                        gain_corrected_results[:, 1] = gain_corrected_results[:, 1]/factor
                        gain_corrected_results[:, 2] = gain_corrected_results[:, 2]/factor
                        gain_corrected_results[:, 3] = gain_corrected_results[:, 3]/factor
                        output_file.create_dataset("gain_corrected", data=gain_corrected_results)

                    if "amplitude_corrected_not_smooht" in output_file and "gain_corrected_smoothed" not in output_file:
                        gain_corrected_smoothed_results = output_file["amplitude_corrected_not_smooht"][()]
                        gain_corrected_smoothed_results[:, 1] = gain_corrected_smoothed_results[:, 1] / factor
                        gain_corrected_smoothed_results[:, 2] = gain_corrected_smoothed_results[:, 2] / factor
                        gain_corrected_smoothed_results[:, 3] = gain_corrected_smoothed_results[:, 3] / factor
                        output_file.create_dataset("gain_corrected_smoothed", data=gain_corrected_smoothed_results)

                    output_file.close()
                    #break

                    '''
                    for i in range(0, len(source_velocities)):
                        results_data[experiment]['polarizationU1'][i][1] = \
                            results_data[experiments]['polarizationU1'][i][1] / factor

                        results_data[experiment]['polarizationU9'][i][1] = \
                            results_data[experiment]['polarizationU9'][i][1] / factor

                        results_data[experiment]['polarizationAVG'][i][1] = \
                            results_data[experiment]['polarizationAVG'][i][1] / factor
                    '''

        '''
        with open(result_files_dir + result_file, "w") as results_out:
            results_out.write(json.dumps(results_data, indent=2))
        '''
        #break

    plt.show()


if __name__ == "__main__":
    main()
    sys.exit(0)
