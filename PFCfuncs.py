import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


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
