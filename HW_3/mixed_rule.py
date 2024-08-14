from sympy import exp, symbols, lambdify, sin, gamma
import numpy as np
import matplotlib.pyplot as plt


class MixedRule:

    def __init__(self, point_a, point_b, T, points, expected_value):
        self.a = point_a
        self.b = point_b
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
        return abs(self.tN - self.t0) / N

    def apply_mixed_rule(self):
        t = symbols("t")
        x = exp(t - exp(-t))
        dx_dt = (1 + exp(-t)) * x
        f_x_t = self.curr_func(x)
        new_integrand = f_x_t * dx_dt
        print("NEW INTEGRAND: ", new_integrand)
        return new_integrand

    def apply_trapezoidal_rule(self, h, f):
        t = symbols("t")
        result = 0
        start = self.t0
        end = self.tN - h
        f_0 = f.evalf(subs={t: end})
        f_N_1 = f.evalf(subs={t: self.tN})
        while start < end - h:
            f_i = f.evalf(subs={t: start + h})
            result += f_i
            start += h
        result = h * result + (h / 2) * (f_0 + f_N_1)
        return result

    def get_simpson_estimation(self):
        f = self.apply_mixed_rule()
        h_N = self.get_optimal_h(self.N)
        S_N = self.apply_trapezoidal_rule(h_N, f)
        h_2N = self.get_optimal_h(2 * self.N)
        S_2N = self.apply_trapezoidal_rule(h_2N, f)
        S = ((4 / 3) * S_2N) - ((1 / 3) * S_N)
        print("EXPECTED: ", self.expected_value)
        print("FINAL RESULT: ", S)
        return S

    def plot_function(
        self,
        function=None,
        symbol=None,
        start=None,
        end=None,
        values=None,
        new_interval=False,
    ):
        if function and values:
            t_range = np.linspace(start, end)
            x = symbols("x")
            f = lambdify(x, function(x))
            x_values, y_values = values
            plt.yscale("log")
            plt.plot(x_values, y_values, "o", label="Desviación")
            plt.plot(t_range, f(t_range), "-", label="1/N^4")
            plt.xlabel("N")
            plt.ylabel("Desviación(N)")
            plt.title("Desviación(N)")

            plt.title(f"{self.curr_func.__name__}")
        elif values:
            x_values, y_values = values
            plt.yscale("log")
            plt.plot(x_values, y_values, "o", label="Desviación")
            plt.xlabel("N")
            plt.ylabel("Desviación(N)")
            plt.title("Desviación(N)")

        elif function:
            print("ENTRA AQUI")
            t_range = np.linspace(start, end)
            if new_interval:
                t = symbols("t")
                coord = [function.subs(t, i) for i in t_range]
                plt.plot(t_range, coord, "-")
                plt.title(f"{function}")
            else:
                x = symbols("x")
                f = lambdify(x, function(x))
                plt.plot(t_range, f(t_range), "-")
                plt.title(f"{function.__name__}")
            plt.xlabel(symbol)
            plt.ylabel(f"f({symbol})")

        plt.legend()
        plt.show()

    def get_error_deviation(self):
        N_values = [i for i in range(1, self.N)]
        s_values = []
        error_values = []
        expected_value = self.expected_value
        for N in N_values:
            self.N = N
            s = self.get_simpson_estimation()
            s_values.append(s)
            deviation = abs(s - expected_value)
            error_values.append(deviation)
        return N_values, error_values


def func_4_5_17(x):
    return (x ** (-3 / 2)) * sin(x / 2) * exp(-x)


# result = MixedRule(
#     point_a=0,
#     point_b=np.inf,
#     T=func_4_5_17,
#     points=30,
#     expected_value=(np.pi * ((5 ** (1 / 2)) - 2)) ** (1 / 2),
# )
# result.plot_function(function=func_4_5_17, symbol="x", start=1, end=100)

# new_integrand = result.apply_mixed_rule()
# result.plot_function(
#     function=new_integrand,
#     symbol="t",
#     start=-4,
#     end=4,
#     values=None,
#     new_interval=True,
# )
# result.t0 = -4
# result.tN = 4
# N_values, error_values = result.get_error_deviation()
# print(N_values, error_values)
# result.plot_function(
#     function=lambda x: 1 / x**4,
#     symbol="N",
#     start=1,
#     end=len(N_values),
#     values=(N_values, error_values),
# )


def func_4_5_18(x):
    return x ** (-2 / 7) * exp(-(x**2))


result = MixedRule(
    point_a=0,
    point_b=np.inf,
    T=func_4_5_18,
    points=30,
    expected_value=(1 / 2) * gamma(5 / 14),
)
# result.plot_function(function=func_4_5_18, symbol="x", start=0, end=100)

# new_integrand = result.apply_mixed_rule()
# result.plot_function(
#     function=new_integrand,
#     symbol="t",
#     start=-4,
#     end=3,
#     values=None,
#     new_interval=True,
# )
result.t0 = -4
result.tN = 3
N_values, error_values = result.get_error_deviation()
print(N_values, error_values)
result.plot_function(
    function=lambda x: 1 / x**4,
    symbol="N",
    start=1,
    end=len(N_values),
    values=(N_values, error_values),
)
