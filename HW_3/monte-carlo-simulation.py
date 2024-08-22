import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
from scipy.constants import k


class IsingModel:

    def __init__(self, dimension):
        self.N = dimension
        self.matrix_neighbor = None

    def get_N_matrix(self, dimension):
        arr = []
        for i in range(dimension):
            row = [
                1 if random.uniform(0, 1) <= 0.7 else -1
                for i in range(dimension)
            ]
            arr.append(row)
        df = pd.DataFrame(arr)
        return df

    def redefine_sites(self, dimension):
        sites = []
        for i in range(1, dimension + 1):
            row = [dimension * (i - 1) + j for j in range(1, dimension + 1)]
            sites.append(row)
        df = pd.DataFrame(sites)
        df = df.sort_index(ascending=False)
        df.index = [i for i in range(dimension)]
        return df

    def get_matrix_neighbor(self, dimension):
        original_sites = self.redefine_sites(dimension)
        row, column = original_sites.shape
        idx = []
        nearest_sites = []
        for i in range(row):
            for j in range(column):
                site = original_sites[j][i]
                up = (
                    site + dimension
                    if site + dimension <= dimension**2
                    else site - dimension * (dimension - 1)
                )
                up_row, up_column = np.where(original_sites == up)
                right = (
                    site + 1 if site % dimension != 0 else site - dimension + 1
                )
                right_row, right_column = np.where(original_sites == right)
                down = (
                    site - dimension
                    if site - dimension >= 1
                    else site + dimension * (dimension - 1)
                )
                down_row, down_column = np.where(original_sites == down)
                left = (
                    site - 1
                    if (site - 1) % dimension != 0
                    else site + dimension - 1
                )
                left_row, left_column = np.where(original_sites == left)
                idx.append(f"{(i, j)}")
                nearest_sites.append(
                    [
                        (up_row[0], up_column[0]),
                        (right_row[0], right_column[0]),
                        (down_row[0], down_column[0]),
                        (left_row[0], left_column[0]),
                    ]
                )
        matrix_neighbor = pd.DataFrame(
            nearest_sites, index=idx, columns=["up", "right", "down", "left"]
        )
        self.matrix_neighbor = matrix_neighbor

    def execute_metropolis_algorithm(
        self, original_matrix, matrix_neighbor, row, column
    ):
        J = 0.5
        h = 0
        new_configuration = original_matrix.copy(deep=True)
        new_configuration.loc[row, column] *= -1
        coord = matrix_neighbor.loc[f"{(row, column)}"]
        up_row, up_column = coord["up"]
        up_value = new_configuration.loc[(up_row, up_column)]
        right_row, right_column = coord["right"]
        right_value = new_configuration.loc[(right_row, right_column)]
        down_row, down_column = coord["down"]
        down_value = new_configuration.loc[(down_row, down_column)]
        left_row, left_column = coord["left"]
        left_value = new_configuration.loc[(left_row, left_column)]
        delta_energy = (
            2
            * J
            * new_configuration.loc[(row, column)]
            * (up_value + right_value + down_value + left_value)
            + 2 * h * new_configuration.loc[(row, column)]
        )
        T = 3e23  # This constant should be dynamic
        probability = min(np.exp(-delta_energy / k * T), 1)
        acceptance = True
        if probability != 1:
            r = random.uniform(0, 1)
            if r > probability:
                acceptance = False
        if acceptance:
            return new_configuration, acceptance
        return None, acceptance

    def define_configurations(self, initial_configuration, sweeps):

        configurations = [
            initial_configuration,
        ]
        for i in range(sweeps):
            for i in range(self.N):
                for j in range(self.N):
                    last_config = configurations[-1]
                    candidate, acceptance = self.execute_metropolis_algorithm(
                        last_config, self.matrix_neighbor, i, j
                    )
                    if acceptance:
                        configurations.append(candidate)
            print("CURR CONFiG: ", len(configurations))
        return configurations

    def get_montecarlo_internal_energy(self, initial_configuration, sweeps):
        configurations = self.define_configurations(
            initial_configuration, sweeps
        )
        J = 0.5
        h = 0

        measurements = []
        sum_energy = 0
        energies = []
        counter = 0
        for configuration in configurations:
            rows, columns = configuration.shape
            for i in range(rows):
                for j in range(columns):
                    coord = self.matrix_neighbor.loc[f"{(i, j)}"]
                    up_row, up_column = coord["up"]
                    up_value = configuration.loc[(up_row, up_column)]
                    right_row, right_column = coord["right"]
                    right_value = configuration.loc[(right_row, right_column)]
                    down_row, down_column = coord["down"]
                    down_value = configuration.loc[(down_row, down_column)]
                    left_row, left_column = coord["left"]
                    left_value = configuration.loc[(left_row, left_column)]
                    energy = (
                        -J
                        * configuration.loc[(i, j)]
                        * (up_value + right_value + down_value + left_value)
                        - h * configuration.loc[(i, j)]
                    )
                    sum_energy += energy
            counter += 1
            measurements.append(counter)
            print("SUM_ENERGY: ", (sum_energy / (len(measurements))))
            energies.append((sum_energy / (len(measurements))))
        return (energies, measurements)  ## ERROR CALCULATION IS MISSING.

    def define_sweeps(self):
        initial_configuration = self.get_N_matrix(16)
        energy_arr = []
        measurements_arr = []
        sweeps = 100
        energies, measurements = self.get_montecarlo_internal_energy(
            initial_configuration, sweeps
        )
        energy_arr.extend(energies)
        measurements_arr.extend(measurements)
        return measurements_arr, energy_arr

    def plot_expected_value(self):
        measurements_arr, energy_arr = self.define_sweeps()
        t = np.linspace(0, len(measurements_arr))
        plt.plot(measurements_arr, energy_arr)
        # plt.plot(t, 0 * t)
        plt.xlabel("measurements")
        plt.ylabel("<E>")
        plt.legend()
        plt.show()


result = IsingModel(16)
result.get_matrix_neighbor(16)
expected_energy = result.plot_expected_value()
