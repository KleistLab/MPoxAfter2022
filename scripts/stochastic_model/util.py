import numpy as np

def ab_concentration(t, x0, f, delta_s, delta_l):

    xs = x0 * f * np.exp(-delta_s * t)
    xl = x0 * (1-f) * np.exp(-delta_l * t)

    return(xs+xl)


def efficacy(c_t, ic50):
    return(c_t/(c_t + ic50)) 