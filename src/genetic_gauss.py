import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.modeling.fitting import LevMarLSQFitter
from functools import reduce
from multiprocessing import Pool
import random
import pandas as pd
import json


test_data_file = "/mnt/WORK/maser/DataProcessingForMaserObservation/output/NotSmooht/cepa/6668/cepa_00_01_03_15_Nov_2019_IRBENE_418.dat"
velocity = np.loadtxt(test_data_file, usecols=(0,), unpack=True)
y_data = np.loadtxt(test_data_file, usecols=(3,), unpack=True)


def compute_gauss(population):
    gauss_lines = population["velocities"]
    amplitudes = population["amplitudes"]
    indexes = [(np.abs(velocity - float(line))).argmin() for line in gauss_lines]
    gaussian = [models.Gaussian1D(amplitudes[index], gauss_lines[index], 0.05, bounds={'stddev': (None, 0.15)}) for index in range(0, len(amplitudes))]
    gg_init = reduce(lambda a, b: a+b, gaussian)
    fitting.SLSQPLSQFitter()
    fit = LevMarLSQFitter()
    gg_fit = fit(gg_init, velocity, y_data)
    return gg_fit


def generate_initial_populations(population_size):
    populations = [] 
   
    for i in range(0, population_size):
        population = {"velocities": [], "amplitudes": []}
        velocity_tmp = []
        cepa_velocity = [-1.77, -2.41, -3.66, -4.01, -4.67]
        tmp = np.random.random()

        if 0.1 < tmp < 0.2:
            initial_population_size = len(cepa_velocity)
        elif 0.3 < tmp < 0.4:
            initial_population_size = len(cepa_velocity) + 1
        elif 0.4 < tmp < 0.5:
            initial_population_size = len(cepa_velocity) + 2
        elif 0.5 < tmp < 0.6:
            initial_population_size = len(cepa_velocity) + 3
        elif 0.6 < tmp < 0.7:
            initial_population_size = len(cepa_velocity) + 4
        elif 0.7 < tmp < 0.8:
            initial_population_size = len(cepa_velocity) + 5
        elif 0.8 < tmp < 0.9:
            initial_population_size = len(cepa_velocity) + 6
        elif 0.9 < tmp < 1.0:
            initial_population_size = len(cepa_velocity) + 7
        else:
            initial_population_size = len(cepa_velocity)

        if initial_population_size == len(cepa_velocity):
            velocity_tmp = cepa_velocity

        else:
            velocity_tmp = cepa_velocity
            for p in range(0, initial_population_size):
                tmp_index = random.randint(0, len(velocity_tmp) - 1)
                velocity_tmp.insert(tmp_index, velocity_tmp[tmp_index] + tmp)

        velocity_tmp = sorted(velocity_tmp, reverse=True)
        indexes = [(np.abs(velocity - float(line))).argmin() for line in velocity_tmp]
        amplitudes = [max(y_data[index - 5:index + 5]) for index in indexes]
        population = {"velocities": velocity_tmp, "amplitudes": amplitudes}
        populations.append(population)

    return populations


def fit(gg_fit):
    return np.sqrt(np.sum(np.abs(y_data**2 - gg_fit(velocity)**2)))


def fitness_evaluation(populations):
    p = Pool(12)
    gg_fits = p.map(compute_gauss, populations)
    fitness = p.map(fit, gg_fits)
    p.close()
    p.terminate()
    #p.join()
    return fitness


def select_elite(populations, fitness):
    selected_individuals_count = len(populations)//2
    if selected_individuals_count % 2 != 0:
        selected_individuals_count += 1

    df = pd.DataFrame({'fitness': fitness})
    fitness_of_selected_individuals = df.nsmallest(selected_individuals_count, 'fitness')
    selected_individuals_indexes = list(fitness_of_selected_individuals.to_dict()["fitness"].keys())
    selected_individuals_indexes = sorted(selected_individuals_indexes)

    selected_individuals = [populations[p] for p in selected_individuals_indexes]
    return selected_individuals


def pairing(elite):
    parents = []

    i = 0
    j = 1

    while j < len(elite):
        parents.append([elite[i], elite[j]])
        i += 1
        j += 1


    return parents


def mutations(parents, max_line_count):
    new_generations = []
    cepa_velocity = [-1.77, -2.41, -3.66, -4.01, -4.67]
    min_velocity_count = len(cepa_velocity)
    max_line_count = max_line_count

    for parent in parents:

        if len(parent[0]["velocities"]) == min_velocity_count and len(parent[1]["velocities"]) == min_velocity_count:
            new_generation = {"velocities": parent[0]["velocities"], "amplitudes": parent[0]["amplitudes"]}

        elif len(parent[0]["velocities"]) == min_velocity_count and len(parent[1]) > min_velocity_count:
            new_generation =  {"velocities": parent[1]["velocities"], "amplitudes": parent[1]["amplitudes"]}

        elif len(parent[1]["velocities"]) == min_velocity_count and len(parent[0]) > min_velocity_count:
            new_generation =  {"velocities": parent[0]["velocities"], "amplitudes": parent[0]["amplitudes"]}

        else:
            new_generation = {"velocities": cepa_velocity}

            new_generation["velocities"].extend(parent[0]["velocities"])
            new_generation["velocities"].extend(parent[1]["velocities"])
            
            new_generation["velocities"] = list(set(new_generation["velocities"]))
            indexes = [(np.abs(velocity - float(line))).argmin() for line in new_generation["velocities"]]
            amplitudes = [max(y_data[index - 5:index + 5]) for index in indexes]
            new_generation["amplitudes"] = amplitudes
            
        tmp = np.random.random()
        tmp_index = random.randint(0, len(new_generation["velocities"]) -1)

        if 0.5 < tmp < 0.75:
            new_generation["velocities"][tmp_index] = new_generation["velocities"][tmp_index] + tmp
            new_generation["amplitudes"][tmp_index] = new_generation["amplitudes"][tmp_index] + tmp

        elif 0.75 <= tmp <= 1.0:
            new_generation["velocities"][tmp_index] = new_generation["velocities"][tmp_index] - tmp
            new_generation["amplitudes"][tmp_index] = new_generation["amplitudes"][tmp_index] - tmp
       
        new_generation["velocities"] = sorted(new_generation["velocities"], reverse=True)

        if len(new_generation["velocities"]) > max_line_count:
            new_generation = {"velocities":new_generation["velocities"][0:max_line_count], "amplitudes":new_generation["amplitudes"][0:max_line_count]}
        new_generations.append(new_generation)

    return new_generations


def plot_best_individual(populations, label, text):
    best_gauss_lines = populations["velocities"][0]
    gg_fit = compute_gauss(best_gauss_lines)
    plt.plot(velocity, y_data, "r-", label="original data")
    plt.plot(velocity, gg_fit(velocity), "g*", label="modulate data")
    plt.xlim((-5, -1))
    plt.text(-5, 400, text)
    plt.legend()
    plt.savefig(label)
    plt.close('all')


def run(population_size, max_line_count):
    tmp = 100
    generations = 1000
    population_size = population_size
    results = dict()

    for r in range(0, tmp):

        print("generation", 0)
        populations = generate_initial_populations(population_size) 
        fitness = fitness_evaluation(populations)
       

        selected_elite_ = select_elite(populations, fitness)
        parents = pairing(selected_elite_)

        i = 0
        for gen in range(0, generations):
            print("generation", i + 1)
            if len(populations) < 2:
                break

            if min(fitness) < 0.0000001:
                break

            populations = mutations(parents, max_line_count)
            fitness = fitness_evaluation(populations)
            selected_elite_ = select_elite(populations, fitness)
            parents = pairing(selected_elite_)

            i += 1

        p = [str(pi) for pi in populations[0]]
        text = "_\n".join(p) + "\n_" + str(fitness[0])
        results[r] = {"best_fitness":fitness[0], "velocities":p}
        label = str(r) + ".png"
        # plot_best_individual(populations, label, text)
  
    return results


def main():
    params = [(10, 10), (10, 15), (10, 20), (10, 30), (15, 10), (15, 15), (15, 20), (15, 30)]

    for param in params:
        results = run(param[0], param[1])

        with open(str(param[0]) + "_" + str(param[1]) + ".json", "w") as result_file:
            result_file.write(json.dumps(results, indent=2))


    sys.exit(0)


if __name__ == "__main__":
    main()

