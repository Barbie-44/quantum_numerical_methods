import numpy as np
from sympy import symbols
import matplotlib.pyplot as plt
from HW_1.problem_1.commons import Coordinates


class LeastSquares(Coordinates):

    def get_chevyshev_polynomials_matrix(self, x_coordinates):
        """
        REMARK: This function only considers the first 4 Chevyshev polynomials.
        """
        polynomial_matrix = []
        for x in x_coordinates:
            T_0 = 1
            T_1 = x
            T_2 = (2 * x * T_1) - T_0
            T_3 = (2 * x * T_2) - T_1
            polynomial_matrix.append([T_0, T_1, T_2, T_3])
        return np.array(polynomial_matrix)

    def get_c_coefficients(self):
        x_coordinates, y_coordinates = self.set_coordinates(self.end)
        polynomial_matrix = self.get_chevyshev_polynomials_matrix(
            x_coordinates
        )
        matrix_transpose = polynomial_matrix.transpose()
        matrix_multiplication = np.matmul(matrix_transpose, polynomial_matrix)
        invertible_matrix = np.linalg.inv(matrix_multiplication)
        print("INVERTIBLE_MATRIX: ", invertible_matrix)
        coefficients_array = np.matmul(
            invertible_matrix,
            np.matmul(matrix_transpose, np.array(y_coordinates)),
        )
        return coefficients_array

    def get_final_function(self, x_symbol):
        coefficients = list(self.get_c_coefficients())
        print("COEFFICIENTS: ", coefficients)
        T_0 = 1
        T_1 = x_symbol
        T_2 = (2 * x_symbol * T_1) - T_0
        T_3 = (2 * x_symbol * T_2) - T_1
        final_function = (
            coefficients[0] * T_0
            + coefficients[1] * T_1
            + coefficients[2] * T_2
            + coefficients[3] * T_3
        )
        print("FINAL FUNCTION: ", final_function)
        return final_function

    def plot_final_graph(self):
        x = symbols("x")
        linear_eq = self.get_final_function(x)
        print("LINEAR_EQ: ", linear_eq)
        f_xs = [linear_eq.subs(x, i) for i in self.x_coordinates]
        print("X's: ", self.x_coordinates)
        print("F(x)'s: ", f_xs)
        plt.plot(self.x_coordinates, f_xs)
        plt.show()
