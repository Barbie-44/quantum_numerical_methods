import matplotlib.pyplot as plt
from HW_1.problem_1.commons import Coordinates
import numpy as np


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
                (y_coordinates[i + 1] - y_coordinates[i])
                / x_coordinates[i + 1]
                - x_coordinates[i]
            ) * (middle_x_distance - x_coordinates[i])
            final_x_coordinates.extend([x_coordinates[i], middle_x_distance])
            final_y_coordinates.extend([y_coordinates[i], f_x])
        print("X'S: ", final_x_coordinates)
        print("Y's: ", final_y_coordinates)
        return final_x_coordinates, final_y_coordinates

    def get_linear_spline_arrays(self):
        x_array, y_array = self.get_linear_spline()
        return x_array, y_array

    def convine_curves(self):
        self.number_of_points = 8
        x_arr, y_arr = self.get_linear_spline_arrays()
        plt.plot(x_arr, y_arr, "r", label="8 puntos")
        self.number_of_points = 16
        x_arr, y_arr = self.get_linear_spline_arrays()
        plt.plot(x_arr, y_arr, "b", label="16 puntos")

        self.number_of_points = 32
        x_arr, y_arr = self.get_linear_spline_arrays()
        plt.plot(x_arr, y_arr, "g", label="32 puntos")

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
