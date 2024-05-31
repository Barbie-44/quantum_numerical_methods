import numpy as np
from scipy.linalg import eigh


class EigenSolver:
    pass


class Hamiltonian(EigenSolver):
    """
    GOAL:
        Diagonalize the XX-Z Hamiltonian by ED method:
        1. Construct the Hamiltonian for N sites
        2. Get Eigenvalues and eigenstates
        3. Get the ground state
        4. Analize the magnetization m_z (order parameter)
        5. Check if there's a phase transition as a function of g=  J_x/J_z
    """

    def __init__(self) -> None:
        self.N = None
        self.I_matrix = np.array([[1, 0], [0, 1]])

    def get_sigma(self, i, sigma):
        sigma_i = self.I_matrix
        for j in range(i, self.N + i):
            if j == i:
                sigma_i = np.kron(sigma_i, sigma)
                continue
            sigma_i = np.kron(sigma_i, self.I_matrix)
        return sigma_i

    def get_sigma_i_and_i_plus_1(self, i, sigma):
        sigma_i = self.get_sigma(i, sigma)
        sigma_i_plus_1 = self.get_sigma(i + 1, sigma)
        print(
            "SHAPE: ",
            sigma_i.shape,
        )
        print("ARRAY_i: ", sigma_i)
        print("ARRA_I_PLUS_1: ", sigma_i_plus_1)
        return np.matmul(sigma_i, sigma_i_plus_1)

    def hamiltonian_builder(self):
        J_x = 1  # Initial J_x
        J_z = 1  # Initial J_z
        S_x = np.array([[0, 1], [1, 0]])
        S_z = np.array([[1, 0], [0, -1]])
        sigma_contributions = []
        for i in range(1, self.N):
            sigma_mat_mul = self.get_sigma_i_and_i_plus_1(i, S_x)
            sigma_contributions.append(sigma_mat_mul)
        HXX = J_x * sum(sigma_contributions)
        sigma_N = self.get_sigma(self.N, S_x)
        second_term = J_x * np.matmul(np.identity(sigma_N.shape[0]), sigma_N)
        total_sigma_z = []
        for i in range(1, self.N):
            current_sigma_z = self.get_sigma(i, S_z)
            total_sigma_z.append(current_sigma_z)
        HZ = J_z * sum(total_sigma_z)
        return HXX + second_term + HZ

    def get_eigenvalues_eigenvectors(self):
        HXX_Z = self.hamiltonian_builder()
        print("COUNT_NONZERO: ", np.count_nonzero(HXX_Z - np.diag(np.diagonal(HXX_Z))))
        print("HXX_Z: ", HXX_Z)
        eigenvalues, eigenvectors = eigh(HXX_Z)
        # print("EIGENVALUES: ", eigenvalues)
        # print("EIGENVECTORES: ", eigenvectors)


problem_2 = Hamiltonian()
problem_2.N = 3
problem_2.get_eigenvalues_eigenvectors()
