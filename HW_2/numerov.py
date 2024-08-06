"""
GOAL:
    - Get the first 10 eigenstates and their respective energies for a particular potential
      via the Numerov method.

APPROACH:
    1. Define the parametric values: Em, w, hbar, m.
    2. Define the potential V(x).
    3. Calculate the optimum distance value (d).
    4. Construct the matrices A and B.
    5. Define the Kinetic operator.
    6. Construct the Hamiltonian.
    7. Get the eigenvalues and eigenstates of the Hamiltonian.

REMARKS:
    - The only input arbitrary value will be the Em energy

"""

import numpy as np
import sympy
from scipy.optimize import root
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt


class Numerov:
    def __init__(self) -> None:
        self.hbar = 1
        self.m = 1
        self.w = 1
        self.Em = 20  # Defatult energy, this value might change by direct definition.
        self.N = 200
        self.d = self.hbar / 2 * np.pi * np.sqrt(2 * self.m * self.Em)

    def double_pit_potential(self, x):
        return self.w * ((x**4) - 2 * x**2)

    def equation(self, x):
        """
        The only purpose of this function is to construct the equation to find the roots of xt's
        """
        value = self.double_pit_potential(x)
        print(value)
        return value - self.Em

    def find_nearest_n(self):
        xt = root(self.equation, self.Em).x[0]
        self.N = round(2 * (xt / (self.d + 4 * np.pi)))

    def id_matrix(self, n, diag):
        ones = [1 + (0 * i) for i in range(n - abs(diag))]
        return np.diag(ones, diag)

    def define_KE_operator(self):
        # self.N = self.find_nearest_n()
        A = (
            self.id_matrix(self.N, -1)
            - 2 * self.id_matrix(self.N, 0)
            + self.id_matrix(self.N, 1)
        ) / self.d
        B = (
            self.id_matrix(self.N, -1)
            + 10 * self.id_matrix(self.N, 0)
            + self.id_matrix(self.N, 1)
        ) / 12
        return -((B.T) * A) / 2

    def get_Hamiltonian(self):
        KE = self.define_KE_operator()
        x = [
            -self.d * (self.N + 1) / 2 + self.d * i
            for i in range(1, self.N + 1)
        ]
        return KE + np.diag(x)

    def get_eig(self):
        H = self.get_Hamiltonian()
        print(H)
        print(H.shape)
        eigenvalues, eigenstates = eigsh(
            H, k=10
        )  # The is the last step of the original problem
        print("eigenvalues: ")
        print(eigenvalues)
        print("eigenstates: ")
        print(eigenstates)
        print(eigenstates.shape)
        x_min, x_max = -20.0, 20.0

        x = np.linspace(x_min, x_max, eigenstates.shape[0])
        for psi in eigenstates:
            plt.plot(x, psi)
        plt.xlabel("x")
        plt.ylabel("Wavefunction (ψ) + Energy")
        plt.legend()
        plt.title("Eigenstates of the Harmonic Oscillator")
        plt.show()


x = Numerov()
x.get_eig()
