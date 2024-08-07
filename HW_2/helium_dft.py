"""
    GOAL:
        Calculate the total energy of helium atom via a DFT algorithm.
    
    APPROACH:
        - Define a initial Gaussian state.
        - Define the functional for DFT.
        - Calculate Ex energy.
        - Calculate Ec energy.
        - Calculate total energy.
    REMARKS:
        - Calculation of Et, Ev and Ej terms are identical as HF.

    CONTEXT:
        HF algorithm is inefficient in the calculation of the He total energy. DFT algorithm provides a better approach. 
        Total energy is the sum of the kinetic energy (Et), potential energy (Ev), potential energy due to Coulomb repulsion (Ej),
        exchange energy (Ex) and correlation energy (Ec) of electrons.
        """

import numpy as np
from scipy.integrate import quad
from sympy import symbols, diff


class DFT:
    def __init__(self) -> None:
        self.e = 1.60218e-19  # Elementary charge (C)
        self.Z = 1.0  # Nuclear charge
        self.epsilon_0 = 8.85419e-12
        self.ansatz = None

    def ansatz(self):  # This might turn into a lambda method
        """
        This method use the Thomas-Fermi density function for the first iteration of n0.
        """
        if not self.ansatz:
            r = symbols("r")
            a0 = 5.29177210903e-11  # Bohr radius in meters
            rs = (
                (((9 * np.pi**2) / 128) ** (1 / 3)) * a0 * (self.Z * self.e) ** (-1 / 3)
            )  # Screening radius
            self.ansatz = (
                ((2 * (self.Z * self.e)) ** (3 / 2)) / (3 * (np.pi**2) * a0**3)
            ) * ((1 - (r / rs)) ** (3 / 2))

    def solve_v_ext(self, r):
        return -self.Z * self.e**2 / (4 * np.pi * self.epsilon_0 * np.abs(r))

    def solve_v_Hartree(self, n0, r):
        integral, _ = quad(lambda r_prime: n0 / np.abs(r - r_prime), -np.inf, np.inf)
        return integral * self.e**2 / (4 * np.pi * self.epsilon_0)

    def solve_v_xc(self, n0):
        n = symbols("n")
        epsilon_xc = -(3 / 4) * ((3 * np.pi) ** (1 / 3)) * (n ** (1 / 3))
        epsilon_prime = diff(epsilon_xc, n)
        return epsilon_prime.subs(n0)

    def solve_kohn_sham_eq(self):
        pass

    def calculate_he_energy(
        self,
    ):
        n0 = self.anzats()
        if converges:
            pass


x = DFT()
