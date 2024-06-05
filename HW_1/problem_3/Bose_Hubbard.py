import math
import numpy as np
import pandas as pd


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

    def get_arr_dimension(self):
        print(
            "dIMENTION: ",
            (
                math.factorial(self.N_particles + self.M_sites - 1)
                / (math.factorial(self.N_particles) * math.factorial(self.M_sites - 1))
            ),
        )
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
        for i in range(len(v_numbers)):
            pi = np.sqrt(100 * i + 3)
            curr_sum += pi * v_numbers[i]
        return curr_sum

    def assign_functional_tag(self):
        basis_arr = self.get_basis_vectors_array()
        print(basis_arr)
        print("*" * 10)
        counter = 0
        T_tags = {}
        while counter < self.get_arr_dimension():
            curr_sum = 0
            for i in range(1, self.M_sites + 1):
                print("i: ", i)
                pi = np.sqrt(100 * i + 3)
                curr_sum += pi * basis_arr.iloc[counter].loc[f"n{i}"]
                print("basis arra_i: ", basis_arr.iloc[counter].loc[f"n{i}"])
                print("CURRENT SUM: ", curr_sum)

            T_tags[f"v{counter+1}"] = curr_sum
            counter += 1
        T_series = pd.Series(T_tags)
        TSorted = T_series.sort_values(ascending=True)
        print(TSorted)

    def diagonalize_H_kin(self):
        """
        This function diagonalize the kinetic term of Bose-Hubbard Hamiltonian column by column given a set of N v-vectors.
        """
        a_dagger_operator = []
        for i in range(self.M_sites):
            curr_row = [0 for entry in range(self.M_sites)]
            if i != 0:
                curr_row[i - 1] = np.sqrt(i)
            a_dagger_operator.append(curr_row)
        a_dagger_operator = np.array(a_dagger_operator)
        a_operator = a_dagger_operator.T
        print(a_dagger_operator)
        print("TRANSPOSE: ", np.transpose(np.array([2, 1, 0])))
        print("A*v: ", np.matmul(a_operator, np.transpose(np.array([2, 1, 0]))))
        v_components = list(
            np.matmul(a_dagger_operator, np.matmul(a_operator, np.array([2, 1, 0]).T))
        )
        for site in range(1, self.M_sites + 1):
            v_sqrt = np.sqrt(v_components)
        print(self.assign_functional_tag())

        print("v_COMPONENTS: ", v_components)
        v_tag = self.get_T_value(v_components)
        print("V_TAG: ", v_tag)


x = ExactDiagonalization()
x.M_sites = 3
x.N_particles = 3
x.diagonalize_H_kin()
