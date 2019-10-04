import sys
import matplotlib.pyplot as plt
import numpy as np
import multiLinReg as MLR
from scipy import stats
from lmfit import Model
from os import path
import pandas as pd
from polymerModels import WLCmodel, FJCmodel


def outputFiles(dataFiles, addon):
    """Takes a list of file names, removes the rxtension, then adds a new
       extension. returns a list"""
    L = []
    for x in range(len(dataFiles)):
        temp = dataFiles[x]
        L.append(temp[:-4] + addon)
    return L


def smooth(x, window_len, window):
    """smooths an numpy array, using a specified window length and
       windowing function. Returns an ndarray."""
    if x.ndim != 1:
        raise ValueError('smooth only accepts 1-D array')
    if x.size < window_len:
        raise ValueError('input vectoir mut be larger than window size')
    if window_len < 3:
        return x
    if window_len % 2 == 0:
        raise ValueError('Window length must be an odd number.')
    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window must be: 'flat', 'hanning', 'hamming', "
                         "'bartlett', or 'blackman'")

    end = window_len - int(window_len)/2

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]

    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')

    return y[int(window_len-end):int(-window_len+end)]


def findBndPtsExact(points, x_points, int_points):
    """
    """
    bound_pts = np.zeros((5, 2))

    slopeSwitch = True
    bndCnt = 0
    for x in range(points):
        if slopeSwitch:
            if int_points[x] == 0:
                bound_pts[bndCnt, 0] = x_points[x]
                bound_pts[bndCnt, 1] = int_points[x]
                slopeSwitch = False
                bndCnt = bndCnt + 1
        else:
            if abs(int_points[x]) > 0:
                bound_pts[bndCnt, 0] = x_points[x-1]
                bound_pts[bndCnt, 1] = int_points[x-1]
                slopeSwitch = True
                bndCnt = bndCnt + 1
        if bndCnt > 4:
            break

    return bound_pts


def findBndPtsNear(points, x_points, int_points):
    """
    """
    bound_pts = np.zeros((5, 2))

    slopeSwitch = True
    bndCnt = 0
    for x in range(points):
        if slopeSwitch:
            if int_points[x] < 0.2:
                bound_pts[bndCnt, 0] = x_points[x]
                bound_pts[bndCnt, 1] = int_points[x]
                slopeSwitch = False
                bndCnt = bndCnt + 1
        else:
            if abs(int_points[x]) > 0.2:
                bound_pts[bndCnt, 0] = x_points[x-1]
                bound_pts[bndCnt, 1] = int_points[x-1]
                slopeSwitch = True
                bndCnt = bndCnt + 1
        if bndCnt > 4:
            break

    return bound_pts


def returnBoundaries(xdata, ydata, resolution):
    """Determines the bounds of the the different parts of the force spectrum
       using the piezo voltage (ydata) and the time series (xdata).
       Assumes the following characteristics: negative slope retract, zero
       slope hold, positive slope approach curve, zero slope hold, negative
       slope retract, zero slope hold
    """
    numPoints = len(ydata)
    points = resolution

    x_points = np.zeros(points)
    y_points = np.zeros(points)
    int_points = np.zeros(points)
    bound_pts = np.zeros((5, 2))
    VboundsXY = np.zeros((5, 2))
    VboundsI = np.zeros(5, dtype=int)

    for x in range(points):
        low = int(x*numPoints/points)
        high = low + int(numPoints/points)
        if high > numPoints:
            high = numPoints
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            xdata[low:high], ydata[low:high])
        x_points[x] = np.mean(xdata[low:high])
        y_points[x] = slope
        int_points[x] = int(slope)

    bound_pts = findBndPtsExact(points, x_points, int_points)

    if 0 in bound_pts[:, 0]:
        for x in range(points):
            int_points[x] = round(y_points[x], 1)
        bound_pts = findBndPtsNear(points, x_points, int_points)

    bndCnt = 0
    for x in range(numPoints):
        if xdata[x] > bound_pts[bndCnt, 0]:
            VboundsXY[bndCnt, 0] = xdata[x-1]
            VboundsXY[bndCnt, 1] = ydata[x-1]
            VboundsI[bndCnt] = x
            # bndSwitch = False
            bndCnt = bndCnt+1
        if bndCnt > 4:
            break

    return VboundsXY, VboundsI


def plotEverything(skipPLT5, originPt, ruptureI, smooth25, separation,
                   baselineS, baselineI, contactS, contactI, result,
                   timeCh1, distance, retractZ, retractD, VboundsXY, VboundsI):

    x_shift = (baselineI - contactI)/(contactS - baselineS)
    y_shift = contactS * x_shift + contactI

    plt.figure(figsize=(20, 10))
    plt.subplot(2, 3, 1)
    plt.title("Z-position (nm)")
    plt.plot(timeCh1, distance)
    plt.plot(VboundsXY[:, 0], VboundsXY[:, 1], 'ro')
    plt.xlabel("Time (s)")
    plt.ylabel("Z-Position (V)")

    plt.subplot(2, 3, 2)
    plt.title("Fitted Retract")
    plt.plot(retractZ, retractD)
    plt.plot(retractZ, baselineS*retractZ + baselineI)
    plt.plot(retractZ, contactS*retractZ + contactI)
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (nm)")
    plt.axis([min(retractZ)-5, max(retractZ)+5, min(retractD)-10,
              max(retractD)+10])

    plt.subplot(2, 3, 3)
    plt.title("Full Retract")
    plt.plot(retractZ, retractD)
    plt.plot(retractZ - x_shift, retractD - y_shift)
    plt.plot(0, 0, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (nm)")

    plt.subplot(2, 3, 4)
    plt.title("Retract")
    plt.plot(retractZ - x_shift, retractD - y_shift)
    plt.plot(0, 0, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (nm)")
    plt.axis([-100, 10, min(retractD - y_shift)-5, 30])

    if skipPLT5:
        plt.subplot(2, 3, 5)
        plt.title("Fit")
        plt.plot(separation[originPt:ruptureI],
                 smooth25[originPt:ruptureI], 'b.')
        # plt.plot(separation[originPt:ruptureI], result.init_fit, 'k--')
        plt.plot(separation[originPt:ruptureI], result.best_fit, 'r-')
        plt.ylabel("Force (nN)")
        plt.xlabel("Separation (nm)")
    else:
        plt.subplot(2, 3, 5)
        plt.title("Fit")
        plt.plot(separation[originPt:ruptureI],
                 smooth25[originPt:ruptureI], 'b.')
        plt.ylabel("Force (nN)")
        plt.xlabel("Separation (nm)")

    # if skipPLT6:
    #     plt.subplot(2, 3, 6)
    #     plt.title("Fit")
    #     plt.plot(separation[originPt:ruptureI],
    #              smooth25[originPt:ruptureI], 'b.')
    #     # plt.plot(separation[originPt:ruptureI], FJCresult.init_fit, 'k--')
    #     plt.plot(FJCresult.best_fit, smooth25[originPt:ruptureI], 'r-')
    #     plt.ylabel("Force (nN)")
    #     plt.xlabel("Separation (nm)")
    # else:
    #     plt.subplot(2, 3, 6)
    #     plt.title("Fit")
    #     plt.plot(separation[originPt:ruptureI],
    #              smooth25[originPt:ruptureI], 'b.')
    #     plt.ylabel("Force (nN)")
    #     plt.xlabel("Separation (nm)")


def mainAnalysis(x1, srcDir, dstDir, csvDir,
                 dataFiles, dataImg, csvOutput, csvRupture):
    currentfile = dataFiles[x1]
    currentpic = dataImg[x1]
    outputfile = csvOutput[x1]
    ruptureFile = csvRupture[x1]

    dataFile = np.genfromtxt(path.join(srcDir, currentfile), skip_header=1)

    # units: s, V, s, m
    (timeCh1, distance, timeDefl, deflection) = dataFile.T
    deflection = deflection*1000000000  # convert deflection to nm

    # Stiffness
    k_L = 0.034  # N/m

    # Find boundaries of setup, approach, retract regions
    # using z-piezo position
    VboundsXY, VboundsI = returnBoundaries(timeCh1, distance, 500)

    # Rescaled vectors to simplify following functions
    approachT = timeCh1[VboundsI[1]:VboundsI[2]]
    approachZ = 4.77*13*distance[VboundsI[1]:VboundsI[2]]
    approachD = deflection[VboundsI[1]:VboundsI[2]]
    approachF = k_L * approachD
    retractT = timeCh1[VboundsI[3]:VboundsI[4]]
    retractZ = 4.77*13*distance[VboundsI[3]:VboundsI[4]]
    retractD = deflection[VboundsI[3]:VboundsI[4]]
    retractF = k_L * retractD

    # Remove data if deflection out of range
    maxDefl = max(retractD)
    x = 0
    while retractD[x] == maxDefl:
        if x > len(retractD)-1:
            break
        x += 1
    retractT = retractT[x:]
    retractZ = retractZ[x:]
    retractD = retractD[x:]
    retractF = retractF[x:]

    # Fit retract curve to get baseline and contact line
    try:
        contactS, contactI, baselineS, baselineI = MLR.multiLinReg(retractZ,
                                                                   retractD)
    except Exception:
        sys.exc_clear()
        print("File %s failed") % (currentfile)
        plt.plot(retractZ, retractD)
        plt.xlabel("Z-position (nm)")
        plt.ylabel("Deflection (nm)")
        plt.savefig(path.join(dstDir, currentpic))
        plt.close()
        # continue

    x_shift = (baselineI - contactI)/(contactS - baselineS)
    y_shift = contactS * x_shift + contactI

    # Linear Regression on approach/retract regions
    # __1 = slope ; __2 = intercept ; __3 = r_value ;
    # __4 = p_value ; __5 = std_error
    # Setup Speed
    setup1, setup2, setup3, setup4, setup5 = stats.linregress(
        timeCh1[0:VboundsI[0]], distance[0:VboundsI[0]])
    # print "setup v =", abs(setup1*13*4.77), "nm/s"

    # Approach Speed
    appr1, appr2, appr3, appr4, appr5 = stats.linregress(approachT, approachZ)
    # print "approch v =", abs(appr1*13*4.77), "nm/s"

    # Retract Speed
    retr1, retr2, retr3, retr4, retr5 = stats.linregress(retractT, retractZ)
    # print "retract v =", abs(retr1*13*4.77), "nm/s"

    smooth11 = smooth((k_L*(retractD - y_shift)/contactS), 11, 'hanning')
    smooth25 = smooth((k_L*(retractD - y_shift)/contactS), 25, 'hanning')
    smooth55 = smooth((k_L*(retractD - y_shift)/contactS), 55, 'hanning')
    smooth75 = smooth((k_L*(retractD - y_shift)/contactS), 75, 'hanning')

    # Find Rupture Force
    ruptureI = np.argmin(retractD)
    ruptureF = k_L*(retractD[ruptureI] - y_shift)/contactS
    ruptureL = (retractZ[ruptureI] - (retractD[ruptureI] - y_shift) - x_shift)

    for x in range(len(retractZ)):
        if (retractZ[x] - x_shift) < 0:
            originPt = x
            break

    # Fit WLC model to rupture
    separation = (retractZ - (retractD - y_shift) - x_shift)

    skipPLT5 = True
    gmod = Model(WLCmodel)
    gmod.set_param_hint('L_C', value=-60.0)
    gmod.set_param_hint('L_P', value=-0.38, min=-0.42, max=-0.34)
    gmod.set_param_hint('a', value=0.0, min=-10.0, max=10.0)
    gmod.set_param_hint('b', value=0.0, min=-10.0, max=10.0)
    params = gmod.make_params()
    try:
        result = gmod.fit(smooth25[originPt:ruptureI],
                          x=separation[originPt:ruptureI])  # method='cobyla'
    except Exception:
        skipPLT5 = False
        # sys.exc_clear() - no longer in python3
    if skipPLT5:
        x_off = result.params['a'].value
        y_off = result.params['b'].value
        WLC_P = result.params['L_P'].value
        WLC_L0 = result.params['L_C'].value
    else:
        x_off = 0.0
        y_off = 0.0
        WLC_P = 0.0
        WLC_L0 = 0.0

    # # Fit FJC model to rupture
    # skipPLT6 = True
    # FJCmod = Model(FJCmodel)
    # FJCmod.set_param_hint('L0')  # , value = -56.0)
    # FJCmod.set_param_hint('b')  # , value = -3.8, min=-4.0, max=-3.6)
    # # FJCmod.set_param_hint('a', value=0.0, min=-5.0, max=5.0)
    # FJCparams = FJCmod.make_params()
    # try:
    #     FJCresult = FJCmod.fit(separation[originPt:ruptureI],
    #                            x=smooth25[originPt:ruptureI])  # method='cobyla'
    # except Exception:
    #     print("FJC failed")
    #     skipPLT6 = False
    #     sys.exc_clear()
    # if skipPLT6:
    #     # x_off = result.params['a'].value
    #     FJC_L0 = FJCresult.params['L0'].value
    #     FJC_b = FJCresult.params['b'].value
    # else:
    #     # x_off = 0.0
    #     FJC_L0 = 0.0
    #     FJC_b = 0.0

    # Add data to pandas DataFrame
    # df = pd.read_excel(path.join(csvDir, 'dataframe.xlsx'))
    # df.append([currentfile, 1000.0*abs(ruptureF), ruptureL, abs(retr1),
    #            WLC_P, WLC_L0, x_off])
    # df.to_excel(path.join(csvDir, 'dataframe.xlsx'), sheet_name='Sheet1')

    # Output Calculations
    output = np.column_stack((retractZ - x_shift, separation, retractD,
                              k_L*(retractD - y_shift)/contactS, smooth11,
                              smooth25, smooth55, smooth75))
    ruptureOut = np.column_stack((separation[originPt:ruptureI],
                                  smooth25[originPt:ruptureI]))

    csvheader = ("z-position(nm),separation(nm),retractD,Force(nN),Force_11"
                 "(nN),Force_25(nN),Force_55(nN),Force_75(nN),v=%d nm/s"
                 % (abs(retr1)))

    ruptureH = "separation(nm),Force_25(nN)"

    np.savetxt(path.join(csvDir, outputfile), output, header=csvheader,
               comments="", delimiter=',')
    np.savetxt(path.join(csvDir, ruptureFile), ruptureOut, header=ruptureH,
               comments="", delimiter=',')

    # Figures
    plotEverything(skipPLT5, originPt, ruptureI, smooth25, separation,
                   baselineS, baselineI, contactS, contactI, result,
                   timeCh1, distance, retractZ, retractD, VboundsXY, VboundsI)

    plt.savefig(path.join(dstDir, currentpic))
    plt.close()

    print("Completed ", x1+1, " of ", len(dataFiles), " files.")
