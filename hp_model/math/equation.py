import numpy as np


class EquationSolver:
    @classmethod
    def solve_equation(cls, pow_coeffs):
        assert len(pow_coeffs) <= 4
        if len(pow_coeffs) == 3:
            return cls.__solve_quadratic(pow_coeffs)
        elif len(pow_coeffs) == 4:
            return cls.__solve_cubic(pow_coeffs)

    @classmethod
    def __solve_cubic(cls, coeffs):
        def convolve(coeffs):
            a, b, c, d = coeffs
            return [1,
                    0,
                    (3 * a * c - b ** 2) / (3 * a ** 2),
                    (2 * b ** 3 - 9 * a * b * c + 27 * a ** 2 * d) / (27 * a ** 3)]

        coeffs = coeffs
        a, b = 1, 0
        if coeffs[1] != 0:
            a, b = coeffs[:2]
            coeffs = convolve(coeffs)
        _, _, p, q = coeffs
        Q = (p / 3) ** 3 + (q / 2) ** 2
        # print("Q = ", Q)
        if Q < 0:
            al, be = (-q / 2 + 1j * np.sqrt(-Q)) ** (1/3), (-q / 2 - 1j * np.sqrt(-Q)) ** (1/3)
        else:
            al, be = (-q / 2 + np.sqrt(Q)) ** (1/3), (-q / 2 - np.sqrt(Q)) ** (1/3)
        roots = [al + be,
                 -(al + be) / 2 + 1j * (al - be) / 2 * np.sqrt(3),
                 -(al + be) / 2 - 1j * (al - be) / 2 * np.sqrt(3)]
        # select only real roots
        roots = [root.real for root in roots if root.imag == 0]
        roots = np.array(roots) - b / (3 * a)
        # select only real roots > 0
        roots = roots[roots > 0]

        return roots

    @classmethod
    def __solve_quadratic(cls, coeffs):
        a, b, c = coeffs
        D = b ** 2 - 4 * a * c
        if D < 0:
            return []
        elif D == 0:
            return [-b / (2 * a)]
        else:
            return [(-b + np.sqrt(D)) / (2 * a),
                    (-b - np.sqrt(D)) / (2 * a)]
