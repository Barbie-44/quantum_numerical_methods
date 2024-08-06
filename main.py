from HW_1.problem_1.polynomial_interpolation import PolynomialInterpolation
from HW_1.problem_1.spline_interpolation import LinearSpline
from HW_1.problem_1.least_squares import LeastSquares
import numpy as np


def main():
    choice = input(
        "Pick one: (a) Interpolation polynomial, (b) spline interpolation, (c) least_squares "
    )
    if choice == "a":
        exercise_a = PolynomialInterpolation()
        exercise_a.function = "b"
        # exercise_a.number_of_points = 32
        exercise_a.end = 1
        exercise_a.convine_curves()
    elif choice == "b":
        exercise_a = LinearSpline()
        exercise_a.function = "b"
        # exercise_a.number_of_points = 32
        exercise_a.end = 1
        exercise_a.convine_curves()
    elif choice == "c":
        exercise_a = LeastSquares()
        exercise_a.function = "b"
        # exercise_a.number_of_points = 32
        exercise_a.end = 1
        exercise_a.regular_spacing = True
        exercise_a.plot_final_graph()


if __name__ == "__main__":
    main()
