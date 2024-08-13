import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, tanh, sinh, cosh, log, exp, sin


class TrapezoidalRule:
    def __init__(self, point_a, point_b):
        self.a = point_a
        self.b = point_b


class DERule(TrapezoidalRule):
    """
    NOTES:
    - log(x) is base 10?
    - Write explicitly the expression for each function.
    """

    def __init__(self, point_a, point_b, T, points, expected_value):
        super().__init__(point_a, point_b)
        self.s = None  # Corresponds to current value of the integral
        self.curr_func = T
        self.N = points
        self.t0 = None
        self.tN = None
        self.expected_value = expected_value

    def get_optimal_h(self, N):
        """
        The error estimation is pending.
        """
        w = np.pi / 2
        return np.log(2 * np.pi * N * w) / N

    def apply_DErule(self):
        """
        Error calculation is pending
        """
        c = 1.0
        t = symbols("t")
        if abs(self.a) == np.inf or abs(self.b) == np.inf:
            """This rule is applied according to table 4.5.14"""
            x = exp(2 * c * sinh(t))
            dx_dt = 2 * c * cosh(t) * exp(2 * c * sinh(t))
        else:
            x = (1 / 2) * (
                self.b + self.a + ((self.b - self.a) * tanh(c * sinh(t)))
            )
            dx_dt = (
                (1 / 2)
                * (self.b - self.a)
                * (1 / (cosh(c * sinh(t))) ** 2)
                * c
                * cosh(t)
            )
        f_x_t = self.curr_func(x)
        # print(x)
        print(f_x_t)
        # print(dx_dt)
        new_integrand = f_x_t * dx_dt
        print("NEW INTEGRAND: ", new_integrand)
        # t_range = np.linspace(-3, 3)
        # It's necessary to check whether or not the t converges to zero in
        # the edges of the interval.
        # coord = [new_integrand.subs(t, i) for i in t_range]
        # plt.plot(t_range, coord)
        # plt.show()
        return new_integrand

    def apply_trapezoidal_rule(self, h, f):
        t = symbols("t")
        result = 0
        start = self.t0
        end = self.tN
        f_0 = f.evalf(subs={t: -3})
        f_N = f.evalf(subs={t: 3})
        while start < end:
            f_i = f.evalf(subs={t: start + h})
            result += f_i
            start += h
        result = h * result + (h / 2) * (f_0 + f_N)
        return result

    def get_error(self, N):
        k = 1
        return np.exp((-k * N) / np.log(N))

    def get_simpson_estimation(self):
        f = self.apply_DErule()
        h_N = self.get_optimal_h(self.N)
        S_N = self.apply_trapezoidal_rule(h_N, f)
        h_2N = self.get_optimal_h(2 * self.N)
        S_2N = self.apply_trapezoidal_rule(h_2N, f)
        S = ((4 / 3) * S_2N) - ((1 / 3) * S_N)
        print("EXPECTED: ", self.expected_value)
        # print("FINAL_RESULT: ", f"{S} ± {error:.20f}")
        # return f"{S} ± {error}"
        print("FINAL RESULT: ", S)
        return S

    def plot_function(self, function):
        t = np.linspace(-10, 50)
        plt.plot(t, function(t))
        plt.show()


def func_4_5_15(x):
    return log(x) * log(1 - x)


# s_values = []
# N_values = []
# error_values = []
# for i in range(1, 25):
#     expected_value = 2 - ((np.pi**2) / 6)
#     result = DERule(
#         point_a=0,
#         point_b=1,
#         T=func_4_5_15,
#         points=i,
#         expected_value=expected_value,
#     )
#     result.t0 = -3
#     result.tN = 3
#     s = result.get_simpson_estimation()
#     s_values.append(s)
#     N_values.append(i)
#     mean_deviation = (1 / len(s_values)) * sum(
#         [si - expected_value for si in s_values]
#     )
#     error_values.append(mean_deviation)


# print("*" * 100)
# print(N_values)
# print(error_values)
# new_n_values = np.linspace(1, 25)
# plt.yscale("log")
# plt.plot(N_values, error_values, "o", label="Desviación")
# plt.plot(
#     new_n_values,
#     (1 / (new_n_values**4)),
#     "-",
#     label="1/N^4",
# )
# plt.xlabel("N")
# plt.ylabel("Desviación(N)")
# plt.legend()
# plt.title("Desviación(N)")
# plt.show()

# result = DERule(
#     point_a=0,
#     point_b=1,
#     T=func_4_5_15,
#     points=25,
#     expected_value=2 - ((np.pi**2) / 6),
# )
# result.t0 = -3
# result.tN = 3
# result.get_simpson_estimation()


def func_4_5_16(x):
    return 1 / (x ** (1 / 2) * (1 + x))


# s_values = []
# N_values = []
# error_values = []
# for i in range(1, 25):
#     expected_value = np.pi
#     result = DERule(
#         point_a=0,
#         point_b=np.inf,
#         T=func_4_5_16,
#         points=i,
#         expected_value=expected_value,
#     )
#     result.t0 = -3
#     result.tN = 3
#     s = result.get_simpson_estimation()
#     s_values.append(s)
#     N_values.append(i)
#     mean_deviation = (1 / len(s_values)) * sum(
#         [abs(si - expected_value) for si in s_values]
#     )
#     error_values.append(mean_deviation)


# print("*" * 100)
# print(N_values)
# print(error_values)
# new_n_values = np.linspace(1, 25)
# plt.yscale("log")
# plt.plot(N_values, error_values, "o", label="Desviación")

# # plt.plot(
# #     new_n_values,
# #     (1 / (new_n_values**4)),
# #     "-",
# #     label="1/N^4",
# # )
# plt.xlabel("N")
# plt.ylabel("Desviación(N)")
# plt.legend()
# plt.title("Desviación(N)")
# plt.show()

# result = DERule(point_a=0, point_b=np.inf, T=func_4_5_16, points=50)
# result.t0 = -3
# result.tN = 3
# result.plot_function(func_4_5_16)

# result.get_simpson_estimation()
# print("EXPECTED: ", np.pi)


def func_4_5_17(x):
    return (1 / x ** (-3 / 2)) * sin(x / 2) * exp(-x)


s_values = []
N_values = []
error_values = []
for i in range(1, 25):
    expected_value = (np.pi * ((5 ** (1 / 2)) - 2)) ** (1 / 2)
    result = DERule(
        point_a=0,
        point_b=np.inf,
        T=func_4_5_17,
        points=i,
        expected_value=expected_value,
    )
    result.t0 = -4.5
    result.tN = 4
    s = result.get_simpson_estimation()
    s_values.append(s)
    N_values.append(i)
    mean_deviation = (1 / len(s_values)) * sum(
        [abs(si - expected_value) for si in s_values]
    )
    error_values.append(mean_deviation)


print("*" * 100)
print(N_values)
print(error_values)
new_n_values = np.linspace(1, 25)
plt.yscale("log")
plt.plot(N_values, error_values, "o", label="Desviación")

# plt.plot(
#     new_n_values,
#     (1 / (new_n_values**4)),
#     "-",
#     label="1/N^4",
# )
plt.xlabel("N")
plt.ylabel("Desviación(N)")
plt.legend()
plt.title("Desviación(N)")
plt.show()


# result = DERule(point_a=0, point_b=np.inf, T=func_4_5_17, points=50)
# # result.t0 = -3
# # result.tN = 3
# # result.plot_function(func_4_5_17)
# result.apply_DErule()
# # print("EXPECTED: ", np.pi)


# t = np.linspace(0, 50)
# plt.yscale("log")
# plt.plot(t, 1 / t**4)
# plt.show()
