import math
import numpy as np
import pandas as pd
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt


class ExactDiagonalization:
    """
    GOAL:
        Get the ground state and ground energy of the Bose-Hubbard Hamiltonian via Exact Diagonalization (ED).

    STEPS:
        1. Construct a vector basis via lexicographic order according to the number of sites N.
        2. Define a tag function T(v) which is one-to-one and sort it.
        3. Write the Bose-Hubbard Hamiltonian DxD parse matrix with this basis for the two components of the Hamiltonian H_int and H_kin.
        4. Apply Lanczos algorithm to the diagonalized Hamiltonian matrix to ge the ground state and ground energy.
        5. Analize different physical quantities of interest.
    """

    def __init__(self):
        self.M_sites = None
        self.N_particles = None
        self.T_tags = None
        self.ind = None
        self.T_sorted = None
        self.J = None
        self.U = 20

    def get_arr_dimension(self):
        return math.factorial(self.N_particles + self.M_sites - 1) / (
            math.factorial(self.N_particles) * math.factorial(self.M_sites - 1)
        )

    def get_basis_vectors_array(self):
        """
        This function constructs an array to store the basis vectors according to the number of sites (M) and particles (N)
        """
        arr_vectors = []
        dimension = self.get_arr_dimension()
        v = 1
        while v < dimension:
            if len(arr_vectors) == 0:
                first_v = [0 for sites in range(self.M_sites)]
                first_v[0] = self.N_particles
                arr_vectors.append(first_v)
                continue
            last_vector = arr_vectors[-1][: self.M_sites]
            last_vector_copy = last_vector.copy()
            last_vector_copy.reverse()
            curr_k = (
                len(last_vector)
                - list(map(bool, last_vector_copy[1 : self.M_sites])).index(True)
                - 1
            )

            curr_vector = []
            for i in range(1, curr_k + 3):
                if i <= curr_k - 1:
                    curr_vector.append(last_vector[i - 1])
                elif i == curr_k:
                    curr_vector.append(last_vector[i - 1] - 1)
                    if curr_k == self.M_sites:
                        continue
                elif i == curr_k + 1:
                    curr_vector.append(self.N_particles - sum(curr_vector))
                elif i == curr_k + 2:
                    curr_vector.append(0)
            arr_vectors.append(curr_vector[: self.M_sites])
            v += 1
        basis_arr = pd.DataFrame(
            arr_vectors, index=([i for i in range(1, int(dimension + 1))])
        ).T
        basis_arr.index = [f"n{i}" for i in range(1, self.M_sites + 1)]
        return basis_arr.T

    def get_T_value(self, v_numbers):
        curr_sum = 0
        for i in range(1, len(v_numbers) + 1):
            pi = np.sqrt(100 * i + 3)
            curr_sum += pi * v_numbers[i - 1]
        return curr_sum

    def assign_functional_tag(self):
        basis_arr = self.get_basis_vectors_array()
        counter = 0
        T_tags = {}
        while counter < self.get_arr_dimension():
            curr_sum = 0
            for i in range(1, self.M_sites + 1):
                pi = np.sqrt(100 * i + 3)
                curr_sum += pi * basis_arr.iloc[counter].loc[f"n{i}"]

            T_tags[counter + 1] = curr_sum
            counter += 1
        self.T_tags = pd.Series(T_tags)
        self.ind = self.T_tags.index
        self.T_sorted = self.T_tags.sort_values(ascending=True)

    def hopping_term(self, state, i, j):
        """
        Compute the hopping term between two basis states.
        """
        state_copy = state.copy()
        value = None
        if i != j and state[j] != 0:
            state_copy.iloc[j] = state.iloc[j] - 1
            state_copy.iloc[i] = state.iloc[i] + 1
            value = (-1) * self.J * (np.sqrt((state.iloc[i] + 1) * state.iloc[j]))
        return value, state_copy

    def get_H_kin(self):
        """
        This function gets the kinetic term of Bose-Hubbard Hamiltonian column by column given a set of N v-vectors.
        """
        basis_states = self.get_basis_vectors_array()
        self.assign_functional_tag()
        D = int(self.get_arr_dimension())
        H_k = np.zeros((D, D))
        for index, state in basis_states.iterrows():
            for i in range(self.M_sites):
                for j in range(self.M_sites):
                    value, v_state = self.hopping_term(state, i, j)
                    if not state.equals(v_state):
                        Tr = self.get_T_value(v_state)
                        basis_v = self.T_tags.where(self.T_tags == Tr).dropna()
                        u_position = basis_v.index.tolist()[0]
                        H_k[u_position - 1][index - 1] = value
        return H_k

    def get_H_int(self):
        """
        This function deduces the interaction term of the Bose-Hubbard Hamiltonian
        """
        basis_states = self.get_basis_vectors_array()
        D = int(self.get_arr_dimension())
        H_int = np.zeros((D, D))
        for index, state in basis_states.iterrows():
            sum_of_part_in_state = 0
            for i in range(self.M_sites):
                ni = state.iloc[i]
                if ni == 0:
                    continue
                sum_of_part_in_state += ni * (ni - 1)
            H_int[index - 1][index - 1] = (self.U / 2) * sum_of_part_in_state
        return H_int

    def get_total_H(self):
        H_kin = self.get_H_kin()
        H_int = self.get_H_int()
        H_total = np.add(H_kin, H_int)
        print(H_total)
        return H_total

    def get_eigenv(self):
        H = self.get_total_H()
        ground_value, ground_state = eigsh(H, k=1)
        ground_value = ground_value[0]
        ground_state = ground_state.reshape(-1)
        print("EIGENVALUE: ", ground_value)
        print("PSI: ", ground_state)
        return ground_value, ground_state

    def get_fluctuation_on_site(self):
        g_value, g_state = self.get_eigenv()
        a = np.diag(np.sqrt(np.arange(1, len(g_state))), 1)
        a_dagger = a.T
        N = np.dot(a_dagger, a)
        expectation_value_N = np.dot(np.conj(g_state), np.dot(N, g_state))
        expectation_value_N2 = np.dot(np.conj(g_state), np.dot(np.dot(N, N), g_state))
        fluctuation_N = np.sqrt(expectation_value_N2 - expectation_value_N**2)
        print(fluctuation_N)
        return fluctuation_N


x = ExactDiagonalization()
x.M_sites = 3
x.N_particles = 3
U_J = []
fluctuations = []
for i in range(1, 10):
    x.J = i
    U_J.append(x.U / x.J)
    fluctuations.append(x.get_fluctuation_on_site())
print("U_J: ", U_J)
print("FLUCTUATIONS: ", fluctuations)
plt.plot(U_J, fluctuations)
plt.show()
