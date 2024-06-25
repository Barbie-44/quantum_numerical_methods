"""
    GOAL:
        Calculate the eigenstates and eigenvalues for the chain spin of S=1/2 and S=1 via the infinite DMRG algorithm

    APPROACH:
        - Define the Hamiltonian model
        - Built operators
        - Define value m as the max eergy to truncate the basis
        - Start growing blocks adding a site per each block
        - Iterate until m value is reached
        - Apply density matrix truncation:
            1.
            2.
            3.
            4.
            5.
            6.
    REMARKS:
        - ...
    BACKGROUND:
        The Heisenberg spin chain consists of a 1D chain of N sites of particles generlly of 1/2 spin.
        The Hamiltonian describes the nearest neighbor spin-spin interaction 
"""

import numpy as np
from pylanczos import PyLanczos


class SpinChain:

    def __init__(self) -> None:
        self.identity = np.array([[1, 0], [0, 1]])
        self.S_z = np.array([[1, 0], [0, -1]])
        self.S_z_left = np.array([[1, 0], [0, -1]])
        self.S_z_right = np.array([[1, 0], [0, -1]])
        self.S_plus = np.array([[0, 1], [0, 0]])
        self.S_plus_left = np.array([[0, 1], [0, 0]])
        self.S_plus_right = np.array([[0, 1], [0, 0]])
        self.S_minus = np.array([[0, 0], [1, 0]])
        self.S_minus_left = np.array([[0, 0], [1, 0]])
        self.S_minus_right = np.array([[0, 0], [1, 0]])
        self.HB = np.array([[0, 0], [0, 0]])
        self.m = None

    def partial_trace(self, rho, dims, keep):
        """
        Perform the partial trace on a density matrix `rho`.

        Parameters:
        rho (ndarray): The density matrix.
        dims (tuple): The dimensions of the subsystems.
        keep (int): The subsystem to keep (0 for the first subsystem, 1 for the second).

        Returns:
        ndarray: The reduced density matrix.
        """
        dim1, dim2 = dims
        if keep == 0:
            # Trace out the second subsystem
            reduced_rho = np.zeros((dim1, dim1), dtype=complex)
            for i in range(dim2):
                for j in range(dim2):
                    reduced_rho += np.trace(
                        rho.reshape([dim1, dim2, dim1, dim2])[:, i, :, j]
                    ) * np.outer([i], [j])
        elif keep == 1:
            # Trace out the first subsystem
            reduced_rho = np.zeros((dim2, dim2), dtype=complex)
            for i in range(dim1):
                for j in range(dim1):
                    reduced_rho += np.trace(
                        rho.reshape([dim1, dim2, dim1, dim2])[i, :, j, :]
                    ) * np.outer([i], [j])
        return reduced_rho

    def construct_block_H(self):
        H_L = (
            np.kron(self.HB, self.identity)
            + np.kron(self.S_z_left, self.S_z)
            + (1 / 2)
            * (
                np.kron(self.S_plus_left, self.S_minus)
                + np.kron(self.S_minus_left, self.S_plus)
            )
        )
        hl_dim = hr_dim = self.HB.shape[0]
        print(hl_dim)
        print(hr_dim)
        H_R = (
            np.kron(self.identity, self.HB)
            + np.kron(self.S_z, self.S_z_right)
            + (1 / 2)
            * (
                np.kron(self.S_plus, self.S_plus_right)
                + np.kron(self.S_minus, self.S_plus_right)
            )
        )
        print(H_R)
        self.S_z_left = np.kron(np.identity(hl_dim), self.S_z)
        self.S_z_right = np.kron(self.S_z, np.identity(hr_dim))
        self.S_plus_left = np.kron(np.identity(hl_dim), self.S_plus)
        self.S_plus_right = np.kron(self.S_plus, np.identity(hr_dim))
        self.S_minus_left = np.kron(np.identity(hl_dim), self.S_minus)
        self.S_minus_right = np.kron(self.S_minus, np.identity(hr_dim))
        H = (
            np.kron(H_L, np.identity(2 * hr_dim))
            + np.kron(np.identity(2 * hl_dim), H_R)
            + np.kron(self.S_z_left, self.S_z_right)
            + (1 / 2) * np.kron(self.S_plus_left, self.S_minus_right)
            + (1 / 2) * np.kron(self.S_minus_left, self.S_plus_right)
        )
        print(H)
        print(H.shape)
        engine = PyLanczos(H, True, 1)  # Find ground state
        lambda_0, ground_state = engine.run()
        print("Eigenvalue: {}".format(lambda_0))
        print("Eigenvector:")
        print(ground_state)
        dims = (hl_dim, hl_dim)
        reducedDM = self.partial_trace(ground_state, dims, keep=0)
        print(reducedDM)
        diag_trace = PyLanczos(reducedDM, True, self.m)
        eigenvalues, eigenvectors = diag_trace.run()
        print("Eigenvalue: {}".format(eigenvalues))
        print("Eigenvector:")
        print(eigenvectors)
        self.rotate_operators(eigenvectors)

    def rotate_operators(self, eigenmatrix):
        print("*" * 10)
        self.HB = np.kron(eigenmatrix.T, np.kron(eigenmatrix.T, eigenmatrix))
        print(self.HB)


x = SpinChain()
x.m = 2
x.construct_block_H()
