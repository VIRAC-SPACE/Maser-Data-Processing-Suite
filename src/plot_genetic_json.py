import sys
import os
import json
from functools import reduce
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.modeling.fitting import LevMarLSQFitter


def get_all_json_file():
    path = "."
    return [file for file in os.listdir(path) if ".json" in file]


def compute_gauss(gauss_lines, velocity, y_data):
    indexes = [(np.abs(velocity - float(line))).argmin() for line in gauss_lines]
    mons = [max(y_data[index - 5:index + 5]) for index in indexes]
    gaussian = [models.Gaussian1D(mons[index], gauss_lines[index], 0.05, bounds={'stddev': (None, 0.15)}) for index in range(0, len(mons))]

    gg_init = reduce(lambda a, b: a+b, gaussian)
    fitting.SLSQPLSQFitter()
    fit = LevMarLSQFitter()
    gg_fit = fit(gg_init, velocity, y_data)

    return gg_fit


def main():
    test_data_file_before = "/home/janis/Documents/maser/output/NotSmooht/cepa/6668/cepa_23_55_31_14_Nov_2019_IRBENE_417.dat"
    test_data_file_after = "/home/janis/Documents/maser/output/NotSmooht/cepa/6668/cepa_00_06_35_15_Nov_2019_IRBENE_419.dat"

    velocity_before = np.loadtxt(test_data_file_before , usecols=(0,), unpack=True)
    y_data_before = np.loadtxt(test_data_file_before , usecols=(3,), unpack=True)

    velocity_after = np.loadtxt(test_data_file_after, usecols=(0,), unpack=True)
    y_data_after = np.loadtxt(test_data_file_after, usecols=(3,), unpack=True)

    for file in get_all_json_file():
        print("file", file)
        fitness = []
        velocities = []
        with open(file) as data:
            genetic_results = json.load(data)
            for iter in genetic_results:
                fitness.append(genetic_results[iter]["best_fitness"])
                velocities.append(genetic_results[iter]["velocities"])

        gauss_lines = velocities[fitness.index(min(fitness))]
        gg_fit_before = compute_gauss(gauss_lines, velocity_before, y_data_before)
        gg_fit_after = compute_gauss(gauss_lines, velocity_after, y_data_after)

        plt.subplot(2, 2, 1)
        plt.plot(velocity_before, y_data_before, "r-", label="original data before")
        plt.plot(velocity_before, gg_fit_before(velocity_before), "g*", label="modulate data before")
        plt.xlim((-5, -1))
        plt.legend()

        plt.subplot(2, 2, 2)
        plt.xlim((-5, -1))
        plt.plot(velocity_after, y_data_after, "r-", label="original data after")
        plt.plot(velocity_after, gg_fit_after(velocity_after), "g*", label="modulate data after")
        plt.xlim((-5, -1))
        plt.legend()

        plt.savefig(file.split(".")[0] + ".png")
        plt.close('all')

    sys.exit(0)


if __name__ == "__main__":
    main()