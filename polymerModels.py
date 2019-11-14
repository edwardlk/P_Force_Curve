import numpy as np


def WLCmodelFull(x, L_P, L_C, a, b):
    """ Worm-Like Chain model
    """
    kBT = 4.114  # pN nm         kT = 4.114 pN nm
    return (kBT/L_P) * (0.25*(1 - (x-a)/L_C)**(-2) - 0.25 + (x-a)/L_C) + b


def WLCmodelNoXY(x, L_P, L_C):
    """ Worm-Like Chain model
    """
    kBT = 4.114  # pN nm         kT = 4.114 pN nm
    return (kBT/L_P) * (0.25*(1 - (x / L_C))**(-2) - 0.25 + (x / L_C))


def WLCmodelImproved(x, L_P, L_C):
    """Worm-like chain model
    """
    kBT = 4.114  # pN nm
    return (kBT/L_P) * (0.25*(1 - (x / L_C))**(-2) - 0.25 + (x / L_C) - 0.8 * (x / L_C)**(2.15))

def FJCmodel(x, L0, b):
    """
    """
    kBT = 0.0406  # nN m
    return L0*(1/np.tanh(b*x/kBT) - kBT/(b*x))
