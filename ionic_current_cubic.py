import numpy as np


def generate_cubic_current(phi_transmem, a, S):
    ionic = S * (1 / a) * phi_transmem * (phi_transmem - a) * (1 - phi_transmem)
    return np.concatenate((ionic, ionic))

def cubic_nondim_gin(phi_resting, a):
    return 666 * 1 / cubic_deriv(phi_resting, a)


def cubic_deriv(phi_transmem, a):
    return (1 / a) * (-3 * phi_transmem ** 2 + 2 * phi_transmem * (1 + a) - a)


def generate_cubic_integral(phi_resting, a):
    return cubic_nondim_gin(phi_resting, a) / 18
