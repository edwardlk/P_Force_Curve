import sys
import matplotlib.pyplot as plt
import numpy as np
import multiLinReg as MLR
from scipy import stats
from lmfit import Model
from os import path
import pandas as pd
from polymerModels import WLCmodelFull, WLCmodelNoXY, FJCmodel
from collections import Counter


def outputFiles(dataFiles, addon):
    """Takes a list of file names, removes the extension, then adds a new
       extension. returns a list"""
    L = []
    for x in range(len(dataFiles)):
        temp = dataFiles[x]
        L.append(temp[:-4] + addon)
    return L


def fitOutputFiles(rupGuess, addon):
    """Takes a dataframe with list of file names in Col1, removes the
       extension, then adds a new extension. can handle multiple events in a
       single file. returns a list"""
    L = []
    temp = rupGuess.iloc[:, 0].value_counts().to_dict()

    for x in temp:
        for y in range(temp[x]):
            L.append(x[:-10] + 'C' + str(y+1) + addon)
    L.sort()
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


def plotEverything(originPt, baselineS, baselineI, contactS, contactI,
                   retractZ_orig, retractD_orig, smooth3, separation,
                   timeCh1, distance, retractZ, retractD, VboundsXY, VboundsI):
    """ Info
    """

    plt.figure(figsize=(20, 10))
    plt.subplot(2, 3, 1)
    plt.title("Z-position (nm)")
    plt.plot(timeCh1, distance, ',')
    plt.plot(VboundsXY[:, 0], VboundsXY[:, 1], 'ro')
    plt.xlabel("Time (s)")
    plt.ylabel("Z-Position (V)")

    plt.subplot(2, 3, 2)
    plt.title("Fitted Retract")
    plt.plot(retractZ_orig, retractD_orig, ',')
    plt.plot(retractZ_orig, baselineS*retractZ_orig + baselineI)
    plt.plot(retractZ_orig, contactS*retractZ_orig + contactI)
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (nm)")
    plt.axis([min(retractZ_orig)-5, max(retractZ_orig)+5,
              min(retractD_orig)-10, max(retractD_orig)+10])
    plt.grid(True, which="both")

    plt.subplot(2, 3, 3)
    plt.title("Full Retract")
    plt.plot(retractZ_orig, retractD_orig, ',')
    plt.plot(retractZ, retractD, ',')
    plt.plot(0, 0, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (nm)")
    plt.grid(True, which="both")

    plt.subplot(2, 3, 4)
    plt.title("Retract")
    plt.plot(retractZ, retractD, ',')
    plt.plot(0, 0, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (nm)")
    plt.axis([-150, 10, min(retractD)-5, 20])
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(10))
    plt.grid(True, which="both")

    plt.subplot(2, 3, 5)
    plt.title("Retract")
    plt.plot(separation, retractD, ',')
    plt.plot(0, 0, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Separation (nm)")
    plt.axis([-150, 10, min(retractD)-5, 20])
    # plt.gca().xaxis.set_major_locator(plt.MultipleLocator(10))
    plt.grid(True, which="both")

    plt.subplot(2, 3, 6)
    plt.title("Retract")
    plt.plot(retractZ, smooth3, ',b')
    plt.plot(separation, smooth3, ',k')
    plt.plot(0, 0, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Separation (nm)")
    plt.axis([-150, 10, min(smooth3)-5, 20])
    # plt.gca().xaxis.set_major_locator(plt.MultipleLocator(10))
    plt.grid(True, which="both")

    # if skipPLT5:
    #     plt.subplot(2, 3, 5)
    #     plt.title("Fit")
    #     plt.plot(separation[originPt:ruptureI],
    #              smooth25[originPt:ruptureI], 'b.')
    #     # plt.plot(separation[originPt:ruptureI], result.init_fit, 'k--')
    #     plt.plot(separation[originPt:ruptureI], result.best_fit, 'r-')
    #     plt.ylabel("Force (nN)")
    #     plt.xlabel("Separation (nm)")
    # else:
    #     plt.subplot(2, 3, 5)
    #     plt.title("Fit")
    #     plt.plot(separation[originPt:ruptureI],
    #              smooth25[originPt:ruptureI], 'b.')
    #     plt.ylabel("Force (nN)")
    #     plt.xlabel("Separation (nm)")

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
    """ Info
    """
    currentfile = dataFiles[x1]
    currentpic = dataImg[x1]
    outputfile = csvOutput[x1]
    # ruptureFile = csvRupture[x1]

    dataFile = np.genfromtxt(path.join(srcDir, currentfile), skip_header=1)

    # units: s, V, s, m
    (timeCh1, distance, timeDefl, deflection) = dataFile.T
    deflection = deflection*1000000000  # convert deflection to nm

    # Stiffness
    k_L = 1.0  # 0.034  # N/m

    # Find boundaries of setup, approach, retract regions
    # using z-piezo position
    VboundsXY, VboundsI = returnBoundaries(timeCh1, distance, 500)

    # Rescaled vectors to simplify following functions
    approachT = timeCh1[VboundsI[1]:VboundsI[2]]
    approachZ = 4.77*13*distance[VboundsI[1]:VboundsI[2]]
    # approachD = deflection[VboundsI[1]:VboundsI[2]]
    retractT = timeCh1[VboundsI[3]:VboundsI[4]]
    retractZ_orig = 4.77*13*distance[VboundsI[3]:VboundsI[4]]
    retractD_orig = deflection[VboundsI[3]:VboundsI[4]]

    # Remove data if deflection out of range
    maxDefl = max(retractD_orig)
    x = 0
    while retractD_orig[x] == maxDefl:
        if x > len(retractD_orig)-1:
            break
        x += 1
    retractT = retractT[x:]
    retractZ_orig = retractZ_orig[x:]
    retractD_orig = retractD_orig[x:]

    # Fit retract curve to get baseline and contact line
    try:
        (contactS, contactI, baselineS,
         baselineI) = MLR.multiLinReg2(retractZ_orig, retractD_orig)
    except Exception:
        sys.exc_clear()
        print("File %s failed") % (currentfile)
        plt.plot(retractZ_orig, retractD_orig)
        plt.xlabel("Z-position (nm)")
        plt.ylabel("Deflection (nm)")
        plt.savefig(path.join(dstDir, currentpic))
        plt.close()
        # continue

    x_shift = (baselineI - contactI)/(contactS - baselineS)
    y_shift = contactS * x_shift + contactI

    retractZ = retractZ_orig - x_shift
    retractD = retractD_orig - y_shift

    originPt = abs(retractZ).argmin()
    retractD = (retractD - baselineS * retractZ) / (contactS - baselineS)

    retractD = (
        retractD - np.average(retractD[np.argmin(abs(retractZ + 150.0)):
                                       np.argmin(abs(retractZ + 200.0))]))

    separation = retractZ - retractD

    # for x in range(len(retractZ)):
    #     if (retractZ[x]) < 0:
    #         originPt = x
    #         break

    # Linear Regression on approach/retract regions
    # __1 = slope ; __2 = intercept ; __3 = r_value ;
    # __4 = p_value ; __5 = std_error
    # Setup Speed
    setup1, setup2, setup3, setup4, setup5 = (
        stats.linregress(timeCh1[0:VboundsI[0]], distance[0:VboundsI[0]]))
    # print "setup v =", abs(setup1*13*4.77), "nm/s"

    # Approach Speed
    appr1, appr2, appr3, appr4, appr5 = stats.linregress(approachT, approachZ)
    # print "approch v =", abs(appr1*13*4.77), "nm/s"

    # Retract Speed
    retr1, retr2, retr3, retr4, retr5 = stats.linregress(retractT, retractZ)
    # print "retract v =", abs(retr1*13*4.77), "nm/s"

    # Force smoothing
    smDict = {
        'smooth1': 11,
        'smooth2': 21,
        'smooth3': 31,
        'smooth4': 41,
    }
    smooth1 = smooth((k_L*retractD), smDict['smooth1'], 'hanning')
    smooth2 = smooth((k_L*retractD), smDict['smooth2'], 'hanning')
    smooth3 = smooth((k_L*retractD), smDict['smooth3'], 'hanning')
    smooth4 = smooth((k_L*retractD), smDict['smooth4'], 'hanning')

    # Find Rupture Force
    # ruptureI = np.argmin(retractD)
    # ruptureF = k_L*retractD[ruptureI]
    # ruptureL = (retractZ[ruptureI] - (retractD[ruptureI]))

    # Fit WLC model to rupture

    # result = 999.0
    # skipPLT5 = False
    # if skipPLT5:
    #     gmod = Model(WLCmodelFull)
    #     gmod.set_param_hint('L_C', value=-60.0)
    #     gmod.set_param_hint('L_P', value=-0.38, min=-0.42, max=-0.34)
    #     gmod.set_param_hint('a', value=0.0, min=-10.0, max=10.0)
    #     gmod.set_param_hint('b', value=0.0, min=-10.0, max=10.0)
    #     params = gmod.make_params()
    #     try:
    #         result = gmod.fit(smooth2[originPt:ruptureI],
    #                           x=separation[originPt:ruptureI])  # method='cobyla'
    #     except Exception:
    #         skipPLT5 = False
    #         # sys.exc_clear() - no longer in python3
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

    # # Fit FJC model to rupture
    # skipPLT6 = False
    # if skipPLT6:
    #     FJCmod = Model(FJCmodel)
    #     FJCmod.set_param_hint('L0')  # , value = -56.0)
    #     FJCmod.set_param_hint('b')  # , value = -3.8, min=-4.0, max=-3.6)
    #     # FJCmod.set_param_hint('a', value=0.0, min=-5.0, max=5.0)
    #     FJCparams = FJCmod.make_params()
    #     try:
    #         FJCresult = FJCmod.fit(separation[originPt:ruptureI],
    #                                x=smooth2[originPt:ruptureI])  # method='cobyla'
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

    # Add data to pandas DataFrame
    # df = pd.read_pickle(path.join(csvDir, "dummy.pkl"))
    # df.loc[x1] = [currentfile, 1000.0*abs(ruptureF), ruptureL, abs(retr1),
    #               WLC_P, WLC_L0, x_off]
    # df.to_pickle(path.join(csvDir, "dummy.pkl"))

    # Output Calculations
    output = np.column_stack((retractZ, separation, retractD, k_L*retractD,
                              smooth1, smooth2, smooth3, smooth4))
    # ruptureOut = np.column_stack((separation[originPt:ruptureI],
    #                               smooth2[originPt:ruptureI]))

    csvheader = ("z-position(nm),separation(nm),retractD,Force(nN),Force_%d"
                 "(nN),Force_%d(nN),Force_%d(nN),Force_%d(nN),v=%d nm/s"
                 % (smDict['smooth1'], smDict['smooth2'], smDict['smooth3'],
                    smDict['smooth4'], abs(retr1)))

    # ruptureH = "separation(nm),Force_25(nN)"

    np.savetxt(path.join(csvDir, outputfile), output, header=csvheader,
               comments="", delimiter=',')
    # np.savetxt(path.join(csvDir, ruptureFile), ruptureOut, header=ruptureH,
    #            comments="", delimiter=',')

    # Figures
    plotEverything(originPt, baselineS, baselineI, contactS, contactI,
                   retractZ_orig, retractD_orig, smooth3, separation,
                   timeCh1, distance, retractZ, retractD, VboundsXY, VboundsI)
    # plt.show()
    plt.savefig(path.join(dstDir, currentpic))
    plt.close()

    print("Completed ", x1+1, " of ", len(dataFiles), " files.")


def fitGuessPlot(dataFile, xDataCol, yDataCol, minGuessID, minGuessRange,
                 fitStartID, minID):
    """
    """
    plt.plot(dataFile[xDataCol], dataFile[yDataCol], ',')
    plt.axis([-150, 10,
              min(dataFile[yDataCol])-5, min(dataFile[yDataCol])+10])
    plt.plot(dataFile[xDataCol].iloc[[minGuessID]],
             dataFile[yDataCol].iloc[[minGuessID]]+5, 'rx')
    plt.plot(dataFile[xDataCol].iloc[[minGuessID+minGuessRange]],
             dataFile[yDataCol].iloc[[minGuessID]]+5, 'r|')
    plt.plot(dataFile[xDataCol].iloc[[minGuessID-minGuessRange]],
             dataFile[yDataCol].iloc[[minGuessID]]+5, 'r|')
    plt.plot(dataFile[xDataCol].iloc[[minID]],
             dataFile[yDataCol].iloc[[minID]], 'gx')
    plt.plot(dataFile[xDataCol].iloc[[fitStartID]],
             dataFile[yDataCol].iloc[[fitStartID]]-2, 'g|')
    plt.plot(dataFile[fitStartID:minID][xDataCol],
             dataFile[fitStartID:minID][yDataCol], ',r')
    plt.ylabel(yDataCol)
    plt.xlabel(xDataCol)
    plt.grid(True, which="both")


def fitAnalysis(x, srcDir, imgDir, csvDir, rupGuess, dataFiles, rupImg,
                rupOutput):
    """ Info
    """
    # Exit if no rupture to analyze
    if rupGuess.at[x, "min_location"] > 0:
        return

    currentfile = rupGuess.at[x, "filename"][:-4] + ".csv"
    minGuess = rupGuess.at[x, "min_location"]
    fitStart = rupGuess.at[x, "fit_start"]
    currentpic = rupImg[x]
    outputfile = rupOutput[x]

    print('x = ' + str(x) + ' ' + currentfile + '  minGuess' + str(minGuess))

    dataFile = pd.read_csv(path.join(srcDir, currentfile))

    # list of y-data columns from files
    # ['z-position(nm)', 'separation(nm)', 'retractD', 'Force(nN)',
    #  'Force_#1(nN)', 'Force_#2(nN)', 'Force_#3(nN)', 'Force_#4(nN)',
    #  'v=309 nm/s']
    yDataColList = dataFile.columns.values.tolist()[2:8]

    yDataCol = yDataColList[5]

    # Find min near minGuess, get row #'s of min and fitStart
    temp = abs(dataFile['z-position(nm)'] - minGuess)
    minGuessID = temp.idxmin()
    minGuessRange = int(7 * len(dataFile['z-position(nm)']) / (
        dataFile['z-position(nm)'].max() - dataFile['z-position(nm)'].min()))
    minID = dataFile[yDataCol][(minGuessID-minGuessRange):
                               (minGuessID+minGuessRange)].idxmin()

    temp = abs(dataFile['z-position(nm)'] - fitStart)
    fitStartID = temp.idxmin()

    # slice datafile
    fitData = dataFile[fitStartID:minID]

    # fit, may need to rescale data

    plt.figure()
    plt.title(currentfile)
    for x in range(5):
        plt.subplot(2, 3, x+1)
        fitGuessPlot(dataFile, 'z-position(nm)', yDataColList[x+1], minGuessID,
                     minGuessRange, fitStartID, minID)
    plt.savefig(path.join(imgDir, currentpic))
    # plt.show()
    plt.close()
