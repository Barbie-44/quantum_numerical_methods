import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, tanh, sinh, cosh, log, exp, sin, lambdify


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
        return abs(self.tN - self.t0) / N

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
        new_integrand = f_x_t * dx_dt
        print("NEW INTEGRAND: ", new_integrand)
        return new_integrand

    def apply_trapezoidal_rule(self, h, f):
        t = symbols("t")
        result = 0
        start = self.t0
        end = self.tN - h
        f_0 = f.evalf(subs={t: end})
        f_N_1 = f.evalf(subs={t: 3})
        while start < end - h:
            f_i = f.evalf(subs={t: start + h})
            result += f_i
            start += h
        result = h * result + (h / 2) * (f_0 + f_N_1)
        return result

    def get_simpson_estimation(self):
        f = self.apply_DErule()
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


def func_4_5_15(x):
    return log(x) * log(1 - x)


# result = DERule(
#     point_a=0,
#     point_b=1,
#     T=func_4_5_15,
#     points=30,
#     expected_value=2 - ((np.pi**2) / 6),
# )
# result.plot_function(function=func_4_5_15, symbol="x", start=0, end=1)
# new_integrand = result.apply_DErule()
# result.plot_function(
#     function=new_integrand,
#     symbol="t",
#     start=-3,
#     end=3,
#     values=None,
#     new_interval=True,
# )
# result.t0 = -3
# result.tN = 3
# N_values, error_values = result.get_error_deviation()
# print(N_values, error_values)
# result.plot_function(
#     function=lambda x: 1 / x**4,
#     symbol="N",
#     start=1,
#     end=len(N_values),
#     values=(N_values, error_values),
# )


def func_4_5_16(x):
    return 1 / (x ** (1 / 2) * (1 + x))


result = DERule(
    point_a=0,
    point_b=np.inf,
    T=func_4_5_16,
    points=30,
    expected_value=np.pi,
)
# result.plot_function(function=func_4_5_16, symbol="x", start=0, end=100)
# new_integrand = result.apply_DErule()
# result.plot_function(
#     function=new_integrand,
#     symbol="t",
#     start=-4,
#     end=4,
#     values=None,
#     new_interval=True,
# )
result.t0 = -4
result.tN = 4
N_values, error_values = result.get_error_deviation()
print(N_values, error_values)
result.plot_function(
    function=lambda x: 1 / x**4,
    symbol="N",
    start=1,
    end=len(N_values),
    values=(N_values, error_values),
)
