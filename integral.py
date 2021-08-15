import numpy as np

C_m = 1.0
"""membrane capacitance, in uF/cm^2"""
g_Na = 120.0
"""Sodium (Na) maximum conductances, in mS/cm^2"""
g_K = 36.0
"""Postassium (K) maximum conductances, in mS/cm^2"""
g_L = 0.3
"""Leak maximum conductances, in mS/cm^2"""
V_Na = 115.0
"""Sodium (Na) Diffusion potentials, in mV"""
V_K = -12.0
"""Postassium (K) Diffusion potentials, in mV"""
V_L = 10.6
"""Leak current Diffusion potentials, in mV"""

def alpha_n(v):
    return 0.1 * ((10 - v) / (np.exp((10 - v) / 10) - 1))


def alpha_m(v):
    return 0.1 * ((25 - v) / (np.exp((25 - v) / 25) - 1))


def d_alpha_n_dv(v):
    return 0.1 * (-1 * (np.exp((10 - v) / 10) - 1) - (10 - v) * -.1 * np.exp((10 - v) / 10)) / (
            (np.exp((10 - v) / 10) - 1) ** 2)


def d_alpha_m_dv(v):
    return 0.1 * (-1 * (np.exp((25 - v) / 25) - 1) - (25 - v) * -.1 * np.exp((25 - v) / 25)) / (
            (np.exp((25 - v) / 25) - 1) ** 2)


def beta_n(v):
    return 0.125 * np.exp(-v / 80)


def beta_m(v):
    return 4 * np.exp(-v / 18)


def d_beta_n_dv(v):
    return (-.125 / 80) * np.exp(-v / 80)


def d_beta_m_dv(v):
    return (-4 / 18) * np.exp(-v / 18)


def n(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def m(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def dn_dv(v):
    return (d_alpha_n_dv(v) * (alpha_n(v) + beta_n(v)) - (d_alpha_n_dv(v) + d_beta_n_dv(v)) * alpha_n(v)) / (
            alpha_n(v) + beta_n(v)) ** 2


def dm_dv(v):
    return (d_alpha_m_dv(v) * (alpha_m(v) + beta_m(v)) - (d_alpha_m_dv(v) + d_beta_m_dv(v)) * alpha_m(v)) / (
            alpha_m(v) + beta_m(v)) ** 2


def dI_dv(v):
    return g_K * (4 * n(v) ** 3 * dn_dv(v) * v + n(v) ** 4 - 4 * n(v) ** 3 * dn_dv(v) * V_K) + g_Na * (
            3 * m(v) ** 2 * dm_dv(v) * v + m(v) ** 3 - 3 * m(v) ** 2 * dm_dv(v) * V_Na) + g_L


def nondimenionsinalize_gin(v):
    R_m = 1 / dI_dv(v)
    """ not sure if the 666 here is right we can play around with it/ask dr lin"""
    gin_nondim = R_m * 666
    return gin_nondim

def rm(v_resting):
    return 1 / dI_dv(v_resting)


"""run this function to get the integral parameter everything else is a helper function"""
def generate_integral(v_resting):
    """v_resting = -70"""
    gin_nondim = nondimenionsinalize_gin(v_resting)
    return gin_nondim /18

def main():
    print(nondimenionsinalize_gin(-70))
    print(generate_integral(-70))


if __name__ == '__main__':
    main()