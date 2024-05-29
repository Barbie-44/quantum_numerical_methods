import numpy as np
import random


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

    def get_normalized_x_coordinates(self, x_coordinates, requires_normalization):
        if requires_normalization:
            self.x_coordinates = [
                (
                    (
                        (2 * (x - min(x_coordinates)))
                        / (max(x_coordinates) - min(x_coordinates))
                    )
                    - 1
                )
                for x in x_coordinates
            ]

    def get_y_coordinates(self, coordinates, function):
        y_coordinates = []
        for x in coordinates:
            print("X: ", x)
            if function == "a":
                f_x = 2 * np.cos(x) + np.sin(2 * x) + np.sqrt(x)
            else:
                f_x = 2 * np.cos(np.pi * x) + np.sin(1 * np.pi * x) + np.sqrt(np.pi * x)
            y_coordinates.append(f_x)
        print(y_coordinates)
        return y_coordinates

    def set_coordinates(self, end, regular_spacing=True, requires_normalization=False):
        if regular_spacing:
            self.x_coordinates = self.set_regular_spacing(self.number_of_points, end)
        else:
            self.x_coordinates = self.set_random_spacing(self.number_of_points, end)
        # self.get_normalized_x_coordinates(self.x_coordinates, requires_normalization)
        print("X_COORDINATES: ", self.x_coordinates)
        return self.x_coordinates, self.get_y_coordinates(
            self.x_coordinates, self.function
        )
