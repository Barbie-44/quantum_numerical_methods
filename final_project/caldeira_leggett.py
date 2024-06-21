import numpy as np
from itertools import product
import pandas as pd


class CaldeiraLeggett:
    """
    GOAL:
        Get eigenvalues and eigenstates of the Caldeira-Legget Hamiltonian via Exact Diagonalization(ED).
    REMARKS:
        .By simplicity all masses and Planck constant are equal to 1.
        .The number of bath oscillators is finite.
    """

    def __init__(self):
        self.N_sys = None
        self.N_env = None
        self.T_tags = None
        self.ind = None
        self.T_sorted = None
        self.arr_dim = None
        self.w_0 = 1

    def get_basis_vectors_array(self):
        """
        This function creates the basis states
        """
        system_range = range(self.N_sys + 1)
        bath_oscillators = [self.N_sys for item in range(self.N_env)]
        reservoir_ranges = [range(max_n + 1) for max_n in bath_oscillators]
        basis = list(product(system_range, *reservoir_ranges))
        basis_arr = pd.DataFrame(basis, index=([i for i in range(len(basis))])).T
        basis_arr.index = [f"n{i}" for i in range(self.N_env + 1)]
        print(basis_arr.T)
        return basis_arr.T

    def assign_functional_tag(self):
        basis_arr = self.get_basis_vectors_array()
        print(basis_arr)
        self.arr_dim = basis_arr.shape[0]
        counter = 0
        T_tags = {}
        while counter < self.arr_dim:
            curr_sum = 0
            for i in range(0, basis_arr.shape[1]):
                pi = np.sqrt(100 * i + 3)
                curr_sum += pi * basis_arr.iloc[counter].loc[f"n{i}"]

            T_tags[counter] = curr_sum
            counter += 1
        self.T_tags = pd.Series(T_tags)
        print(self.T_tags)
        self.ind = self.T_tags.index
        print(self.ind)
        self.T_sorted = self.T_tags.sort_values(ascending=True)
        print(self.T_sorted)

    def get_H_sys(self):
        w_0 = self.w_0
        H_sys_states = [i for i in range(self.N_sys + 1)]
        basis_states = pd.DataFrame(
            H_sys_states, index=([i for i in range(self.N_sys + 1)])
        )
        print(basis_states)
        H_0 = np.zeros((self.N_sys + 1, self.N_sys + 1))
        print(H_0)
        for index, state in basis_states.iterrows():
            value = state.iloc[0]
            if value == 0:
                value = 1 / np.sqrt(2)
            H_0[index][index] = w_0 * value
        print("*" * 10)
        print(H_0)
        print(np.identity((self.N_sys + 1) ** self.N_env))
        env_dimension = (self.N_sys + 1) ** self.N_env
        result = np.kron(H_0, np.identity(env_dimension))
        print(result)
        print(result.shape)

    def get_H_env(self):
        w_i = 1
        system_range = range(self.N_sys + 1)
        bath_oscillators = [self.N_sys for item in range(self.N_env - 1)]
        reservoir_ranges = [range(max_n + 1) for max_n in bath_oscillators]
        basis = list(product(system_range, *reservoir_ranges))
        basis_arr = pd.DataFrame(basis, index=([i for i in range(len(basis))])).T
        basis_arr.index = [f"n{i}" for i in range(self.N_env)]
        basis_arr = basis_arr.T
        print(basis_arr)
        n_max = (self.N_sys + 1) ** self.N_env
        H_env = np.zeros((n_max, n_max))
        print(H_env)
        print(H_env.shape)
        for index, state in basis_arr.iterrows():
            print(index)
            sum_of_part_in_state = 0
            for i in range(basis_arr.shape[1] - 1):
                sum_of_part_in_state += w_i * state.iloc[i]
            H_env[index][index] = sum_of_part_in_state
        print(H_env)
        result = np.kron(np.identity(self.N_sys + 1), H_env)
        print(result)
        print(result.shape)

    # def get_H_int(self):


x = CaldeiraLeggett()
x.N_env = 3
x.N_sys = 1
x.get_H_env()
