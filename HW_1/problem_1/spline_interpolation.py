import matplotlib.pyplot as plt
from HW_1.problem_1.commons import Coordinates


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




# class CubicSpline(Spacing):
#     pass
