import numpy as np


def WLCmodel(x, L_P, L_C, a, b):
    """ Worm-Like Chain model
    """
    kBT = 4.114  # pN nm         kT = 4.114 pN nm
    return (kBT/L_P) * (0.25*(1 - (x-a)/L_C)**(-2) - 0.25 + (x-a)/L_C) + b


def FJCmodel(x, L0, b):
    """
    """
    kBT = 0.0406  # nN m
    return L0*(1/np.tanh(b*x/kBT) - kBT/(b*x))


# # Fit WLC model to rupture
# def WLCfit(retractZ, retractD, smooth25, y_shift, x_shift, originPt, ruptureI):
#     """Need to determine what gmod.fit returns, """
#     separation = (retractZ - (retractD - y_shift) - x_shift)
#
#     skipPLT5 = True
#     gmod = Model(WLCmodel)
#     gmod.set_param_hint('L_C', value=-60.0)
#     gmod.set_param_hint('L_P', value=-0.38, min=-0.48, max=-0.28)
#     gmod.set_param_hint('a', value=0.0, min=-10.0, max=10.0)
#     gmod.set_param_hint('b', value=0.0, min=-10.0, max=10.0)
#     params = gmod.make_params()
#     try:
#         result = gmod.fit(smooth25[originPt:ruptureI],
#                           x=separation[originPt:ruptureI])  # method='cobyla'
#     except Exception:
#         skipPLT5 = False
#         sys.exc_clear()
#     if skipPLT5:
#         x_off = result.params['a'].value
#         y_off = result.params['b'].value
#         WLC_P = result.params['L_P'].value
#         WLC_L0 = result.params['L_C'].value
#     else:
#         x_off = 0.0
#         y_off = 0.0
#         WLC_P = 0.0
#         WLC_L0 = 0.0
#
#
# # Fit FJC model to rupture
# def FJCfit(retractZ, retractD, smooth25, y_shift, x_shift, originPt, ruptureI):
#     """
#     """
#     separation = (retractZ - (retractD - y_shift) - x_shift)
#
#     skipPLT6 = True
#     FJCmod = Model(FJCmodel)
#     # FJCmod.set_param_hint('L0') #, value = -56.0)
#     # FJCmod.set_param_hint('b') #, value = -3.8, min=-4.0, max=-3.6)
#     # FJCmod.set_param_hint('a', value = 0.0, min=-5.0, max=5.0)
#     FJCparams = FJCmod.make_params()
#     try:
#         FJCresult = FJCmod.fit(separation[originPt:ruptureI],
#                                x=smooth25[originPt:ruptureI])  # method='cobyla'
#     except Exception:
#         print("FJC failed")
#         skipPLT6 = False
#         sys.exc_clear()
#     if skipPLT6:
#         # x_off = result.params['a'].value
#         FJC_L0 = FJCresult.params['L0'].value
#         FJC_b = FJCresult.params['b'].value
#     else:
#         # x_off = 0.0
#         FJC_L0 = 0.0
#         FJC_b = 0.0
