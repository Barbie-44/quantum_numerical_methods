import numpy as np
from sympy import symbols, simplify
import matplotlib.pyplot as plt
from HW_1.problem_1.commons import Coordinates


class PolynomialInterpolation(Coordinates):
    """
    Description:
        This class inherits from Spacing class in order to stablish a regular
        or random spacing.
    Pain Points:
        For n=32 points function behavior is strange.
    """

    def set_polinomic_base(self, x, i, x_coordinates):
        l_i = 1
        for j in range(len(x_coordinates)):
            if j != i:
                l_i *= (x - x_coordinates[j]) / (
                    x_coordinates[i] - x_coordinates[j]
                )
        return l_i

    def set_final_polynomial(self, x_symbol):
        x_coordinates, y_coordinates = self.set_coordinates(self.end)
        polynomyal = 0
        for i in range(len(x_coordinates)):
            L_i = self.set_polinomic_base(x_symbol, i, self.x_coordinates)
            polynomyal += y_coordinates[i] * L_i
        print("POLYNOMIAL: ", simplify(polynomyal))
        return simplify(polynomyal)

    def final_polynomial(self):
        x = symbols("x")
        lagrange_polynomial = self.set_final_polynomial(x)
        f_xs = [lagrange_polynomial.subs(x, i) for i in self.x_coordinates]
        print("X's: ", self.x_coordinates)
        print("F(x)'s: ", f_xs)
        return f_xs

    def convine_curves(self):
        self.number_of_points = 8
        y_8 = self.final_polynomial()
        plt.plot(self.x_coordinates, y_8, "r", label="8 puntos")
        self.number_of_points = 16
        y_16 = self.final_polynomial()
        plt.plot(self.x_coordinates, y_16, "b", label="16 puntos")
        self.number_of_points = 20
        y_20 = self.final_polynomial()

        plt.plot(self.x_coordinates, y_20, "g", label="32 puntos")

        x = np.linspace(0, self.end, 400)
        if self.function == "a":
            plt.plot(
                x,
                2 * np.cos(x) + np.sin(2 * x) + np.sqrt(x),
                "y",
                label="función original",
            )
        else:
            plt.plot(
                x,
                2 * np.cos(np.pi * x)
                + np.sin(1 * np.pi * x)
                + np.sqrt(np.pi * x),
                "y",
                label="función original",
            )
        plt.xlabel("x")
        plt.ylabel("f(x)")
        plt.legend()
        plt.show()
