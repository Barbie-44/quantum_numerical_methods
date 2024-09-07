import numpy as np
from sympy import (
    symbols,
    Rational,
    diff,
    Matrix,
    lambdify,
    assoc_legendre,
    simplify,
    hyper,
    S,
)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import eval_legendre


class NewtonRaphson:

    def __init__(self, eta, guess):
        """
        Guess is an iterable [mu, delta]
        """
        self.eta = eta
        self.x0 = guess

    def get_g0_function(self, z):
        m, d = symbols("m d")
        p_1_2 = hyper([S(-1) / 2, S(3) / 2], [1], (1 - z) / 2)

        g0 = self.eta - (((m**2) + (d**2)) ** (1 / 4)) * p_1_2

        return g0

    def get_g1_function(self, z):
        m, d = symbols("m d")
        p_3_2 = hyper([S(-3) / 2, S(5) / 2], [1], (1 - z) / 2)
        g1 = (
            (self.eta * m)
            + ((((m**2) + (d**2)) ** (3 / 4)) * p_3_2)
            + (4 / 3 * np.pi)
        )
        return g1

    def plot_original_functions(self):
        m_values = np.linspace(-100, 100, 100)
        d_values = np.linspace(-100, 100, 100)
        mu, delta = np.meshgrid(m_values, d_values)
        m, d = symbols("m d")
        g0 = self.get_g0_function()
        print("G0: ", g0)
        Z = lambdify((m, d), g0, "numpy")

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(
            mu,
            delta,
            Z(mu, delta),
            cmap="viridis",
            label="g0",
        )
        g1 = self.get_g1_function()
        Z2 = lambdify([m, d], g1)
        ax.plot_surface(
            mu,
            delta,
            Z2(mu, delta),
            cmap="Oranges",
            label="g1",
        )

        ax.set_xlabel("mu")
        ax.set_ylabel("delta")
        ax.set_title("g0, g1")
        ax.legend()

        plt.show()

    def get_jacobian(self, g0, g1):
        m, d = symbols("m d")
        dg0_dm = simplify(diff(g0, m))
        dg0_dd = simplify(diff(g0, d))
        dg1_dm = simplify(diff(g1, m))
        dg1_dd = simplify(diff(g1, d))
        jacobian = Matrix([[dg0_dm, dg0_dd], [dg1_dm, dg1_dd]])
        # inverse = Matrix([[dg1_dd, -dg0_dd], [-dg1_dm, dg0_dm]])
        # det = 1 / ((dg0_dm * dg1_dd) - (dg0_dd * dg1_dm))
        # inverse_jacobian = det * inverse
        return jacobian

    def calculate_x1(self):
        m, d = symbols("m d")

        x0 = Matrix(self.x0)
        error = 1e-2
        max_evaluations = 100
        counter = 0
        z = -m * (((m**2) + (d**2)) ** (-1 / 2))

        g0 = self.get_g0_function(z)
        g1 = self.get_g1_function(z)
        jacobian = self.get_jacobian(g0, g1)
        while (
            abs(float(g0.evalf(subs={m: float(x0[0]), d: float(x0[1])})))
            > error
        ):

            counter += 1
            print("COUNTER: ", counter)
            if counter >= max_evaluations:
                print("NO HAY")
                return 0, 0
            values = {"m": float(x0[0]), "d": float(x0[1])}
            eval_inverse = jacobian.subs(values).inv()
            g00 = g0.evalf(subs=values)
            g11 = g1.evalf(subs=values)
            print("G0: ", g00, " G1: ", g11)
            g_vector = Matrix([g00, g11])
            print("PRODUCT: ", simplify(eval_inverse * g_vector))
            x0 = simplify(x0 - (eval_inverse * g_vector))
            print("x0: ", float(x0[0]), float(x0[1]))
        print("FINAL G0: ", g0.evalf(subs={m: float(x0[0]), d: float(x0[1])}))
        print("FINAL G1: ", g1.evalf(subs={m: float(x0[0]), d: float(x0[1])}))
        return float(x0[0]), float(x0[1])

    def calculate_error(self):
        pass


guess = [1.0, 0.01]
eta_values = np.linspace(-3, 1)
mu_values = []
delta_values = []
counter = 0
for eta in eta_values:
    print("*" * 100)
    print("ETA: ", eta)
    result = NewtonRaphson(eta, guess)
    m, d = result.calculate_x1()
    mu_values.append(m)
    delta_values.append(d)

print(eta_values)
print("*" * 10)
print(mu_values)
print("*" * 10)
print(delta_values)
plt.plot(eta_values, mu_values, "r", label="mu")
plt.plot(eta_values, delta_values, "b", label="delta")
plt.xlabel("eta")
plt.ylabel("mu, delta")
plt.legend()
plt.show()

# m, d = symbols("m d")
# print("G0 EVAL: ", result.get_g0_function().evalf(subs={m: x1_m, d: x1_d}))
# print("G1 EVAL: ", result.get_g1_function().evalf(subs={m: x1_m, d: x1_d}))
