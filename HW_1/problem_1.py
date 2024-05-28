import numpy as np
import random
import math
from sympy import symbols, simplify
import matplotlib.pyplot as plt


class Spacing:

    def set_regular_spacing(self, N: int, end):
        return list(np.linspace(0, end, N))

    def set_random_spacing(self, N: int, end):
        return [random.uniform(0, end) for number in range(N)]


class Coordinates(Spacing):

    def __init__(self):
        self.function = None
        self.number_of_points = None
        self.end = None
        self.x_coordinates = None

    def get_y_coordinates(self, coordinates, function):
        y_coordinates = []
        print(coordinates)
        for x in coordinates:
            if function == "a":
                f_x = 2 * np.cos(x) + np.sin(2 * x) + np.sqrt(x)
            else:
                f_x = 2 * np.cos(np.pi * x) + np.sin(1 * np.pi * x) + np.sqrt(np.pi * x)
            y_coordinates.append(f_x)
        print(y_coordinates)
        return y_coordinates

    def set_coordinates(
        self,
        end,
        regular_spacing=True,
    ):
        if regular_spacing:
            self.x_coordinates = self.set_regular_spacing(self.number_of_points, end)
        else:
            self.x_coordinates = self.set_random_spacing(self.number_of_points, end)
        return self.x_coordinates, self.get_y_coordinates(
            self.x_coordinates, self.function
        )


class PolynomialInterpolation(Coordinates):
    """
    Description:
        This class inherits from Spacing class in order to stablish a regular or random spacing.
    Pain Points:
        For n=32 points function behavior is strange.
    """

    def set_polinomic_base(self, x, i, x_coordinates):
        l_i = 1
        for j in range(len(x_coordinates)):
            if j != i:
                l_i *= (x - x_coordinates[j]) / (x_coordinates[i] - x_coordinates[j])
        return l_i

    def set_final_polynomial(self, x_symbol):
        x_coordinates, y_coordinates = self.set_coordinates(self.end)
        polynomyal = 0
        for i in range(len(x_coordinates)):
            L_i = self.set_polinomic_base(x_symbol, i, self.x_coordinates)
            polynomyal += y_coordinates[i] * L_i
        print("POLYNOMIAL: ", simplify(polynomyal))
        return simplify(polynomyal)

    def plot_final_polynomial(self):
        x = symbols("x")
        lagrange_polynomial = self.set_final_polynomial(x)
        f_xs = [lagrange_polynomial.subs(x, i) for i in self.x_coordinates]
        print("X's: ", self.x_coordinates)
        print("F(x)'s: ", f_xs)
        plt.plot(self.x_coordinates, f_xs)
        plt.show()


class LinearSpline(Coordinates):
    """The goal of this class is to create a polynomial of first degree"""

    def get_linear_spline(self):
        x_coordinates, y_coordinates = self.set_coordinates(self.end)
        final_x_coordinates = []
        final_y_coordinates = []
        for i in range(len(x_coordinates) - 1):
            middle_x_distance = (
                (x_coordinates[i + 1] - x_coordinates[i]) / 2
            ) + x_coordinates[i]
            print("middle_x_distance: ", middle_x_distance)
            f_x = y_coordinates[i] + (
                (y_coordinates[i + 1] - y_coordinates[i]) / x_coordinates[i + 1]
                - x_coordinates[i]
            ) * (middle_x_distance - x_coordinates[i])
            final_x_coordinates.extend([x_coordinates[i], middle_x_distance])
            final_y_coordinates.extend([y_coordinates[i], f_x])
        print("X'S: ", final_x_coordinates)
        print("Y's: ", final_y_coordinates)
        return final_x_coordinates, final_y_coordinates

    def plot_linear_spline(self):
        x_array, y_array = self.get_linear_spline()
        plt.plot(x_array, y_array)
        plt.show()


# exercise_a = LinearSpline()
# exercise_a.function = "a"
# exercise_a.number_of_points = 32
# exercise_a.end = np.pi
# exercise_a.plot_linear_spline()

# class QuadraticSpline(Spacing):
#     pass

# class CubicSpline(Spacing):
#     pass

# class LeastSquares(Spacing): ...

exercise_a = PolynomialInterpolation()
exercise_a.function = "a"
exercise_a.number_of_points = 32
exercise_a.end = np.pi
exercise_a.plot_final_polynomial()
