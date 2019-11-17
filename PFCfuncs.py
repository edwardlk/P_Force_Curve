import sys
import matplotlib.pyplot as plt
import numpy as np
import multiLinReg as MLR
from scipy import stats
from lmfit import Model
from os import path
import pandas as pd
from polymerModels import WLCmodelNoXY, WLCmodelImproved  # , FJCmodel
from PFCplots import plotEverything, plotFits


def outputFiles(dataFiles, addon):
    """Takes a list of file names, removes the extension, then adds a new
       extension.
       Returns a list"""
    L = []
    for x in range(len(dataFiles)):
        temp = dataFiles[x]
        L.append(temp[:-4] + addon)
    return L


def fitOutputFiles(rupGuess, addon):
    """Takes a dataframe with list of file names in Col1, removes the
       extension, then adds a new extension. Can handle multiple events in a
       single file.
       Returns a list"""
    if not isinstance(rupGuess, pd.DataFrame):
        raise ValueError('fitOutputFiles only accepts dataframe')

    L = []
    temp = rupGuess.iloc[:, 0].value_counts().to_dict()

    for x in temp:
        for y in range(temp[x]):
            L.append(x[:-10] + 'C' + str(y+1) + addon)
    L.sort()
    return L


def smooth(x, window_len, window):
    """Smooths an numpy array, using a specified window length and windowing
       function.
       Returns an ndarray."""
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
    for x in range(1, points):
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
    for x in range(1, points):
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
       slope hold, positive slope approach, zero slope hold, negative
       slope retract, zero slope hold.
       Returns 2 ndarrays"""
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
        int_points[x] = round(slope*2)/2

    # plt.figure()
    # plt.plot(x_points, y_points)
    #
    # plt.figure()
    # plt.plot(x_points, int_points)
    # plt.show()
    # plt.close()

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
    k_L = 0.034  # 0.034  # N/m

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
        x += 1
        if x > len(retractD_orig)-1:
            break

    retractT = retractT[x:]
    retractZ_orig = retractZ_orig[x:]
    retractD_orig = retractD_orig[x:]

    # Fit retract curve to get baseline and contact line
    try:
        (contactS, contactI, baselineS,
         baselineI) = MLR.multiLinReg2(retractZ_orig, retractD_orig)
    except Exception:
        print("File {} failed".format(currentfile))
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

    # Figures
    plotEverything(currentpic, abs(retr1), originPt, baselineS, baselineI,
                   contactS, contactI, retractZ_orig, retractD_orig, smooth3,
                   separation, timeCh1, distance, retractZ, retractD,
                   VboundsXY, VboundsI)
    # plt.show()
    plt.savefig(path.join(dstDir, currentpic))
    plt.close()

    print("Completed ", x1+1, " of ", len(dataFiles), " files.")


def fitToWLCmodel(modelName, ydata, xdata):
    """ Info
    """
    L_CG = max(xdata)
    gmodel = Model(modelName, nan_policy='omit')
    gmodel.set_param_hint('L_C', value=1.5*L_CG, min=0.1*L_CG, max=5.0*L_CG)
    gmodel.set_param_hint('L_P', value=3.0, min=0.0, max=100.0)
    params = gmodel.make_params()
    result = gmodel.fit(ydata, params, x=xdata)
    # print('parameter names: {}'.format(gmodel.param_names))
    # print('independent variables: {}'.format(gmodel.independent_vars))
    return result


def fitPrntList(x, y, vel, currentfile, model1, Ymax, Xmax):
    if model1.params['L_P'].stderr is not None:
        L_Perr = model1.params['L_P'].stderr
    else:
        L_Perr = 0.0
    if model1.params['L_C'].stderr is not None:
        L_Cerr = model1.params['L_P'].stderr
    else:
        L_Cerr = 0.0

    return [currentfile, y, vel, Ymax, Xmax, model1.method, model1.nfev,
            model1.ndata, model1.nvarys, model1.chisqr, model1.redchi,
            model1.aic, model1.bic, model1.best_values['L_P'], L_Perr,
            model1.best_values['L_C'], L_Cerr]


def fitAnalysis(x, srcDir, imgDir, csvDir, rupGuess, dataFiles, rupImg,
                rupOutput):
    """ Info
    """
    # Exit if no rupture to analyze
    currentfile = rupGuess.at[x, "filename"][:-4] + ".csv"
    if rupGuess.at[x, "min_location"] > 0:
        print(currentfile + '  Skip')
        return
    minGuess = rupGuess.at[x, "min_location"]
    fitStart = rupGuess.at[x, "fit_start"]
    currentpic = rupImg[x]
    outputfile = rupOutput[x]

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
    adj = int((minID - fitStartID) / 10)

    fitData_flip = -1 * fitData
    fitData_flipMin = fitData_flip[:adj].agg(['mean'])

    fitData_flipXY = fitData_flip.sub(fitData_flipMin.loc['mean'],
                                      axis='columns')

    # Fit Data.
    # WLCmodelNoXY, WLCmodelImproved
    # model1 - No adjustment to force curve
    # model2 - adjust y-offset
    # model3 - adjust x- and y-offset (move beginning of curve to origin)
    fit_df = pd.read_pickle(path.join(csvDir, "dummy.pkl"))

    model1 = fitToWLCmodel(WLCmodelNoXY, fitData_flip[yDataColList[4]],
                           fitData_flip['z-position(nm)'])
    model2 = fitToWLCmodel(WLCmodelNoXY, fitData_flipXY[yDataColList[4]],
                           fitData_flip['z-position(nm)'])
    model3 = fitToWLCmodel(WLCmodelNoXY, fitData_flipXY[yDataColList[4]],
                           fitData_flipXY['z-position(nm)'])
    modelB1 = fitToWLCmodel(WLCmodelImproved, fitData_flip[yDataColList[4]],
                            fitData_flip['z-position(nm)'])
    modelB2 = fitToWLCmodel(WLCmodelImproved, fitData_flipXY[yDataColList[4]],
                            fitData_flip['z-position(nm)'])
    modelB3 = fitToWLCmodel(WLCmodelImproved, fitData_flipXY[yDataColList[4]],
                            fitData_flipXY['z-position(nm)'])
    vel = int(fitData.columns[-1][2:-4])
    fit_df.loc[len(fit_df)] = fitPrntList(x, 1, vel, outputfile, model1,
                                          fitData_flip[yDataColList[4]].max(),
                                          fitData_flip['z-position(nm)'].max())
    fit_df.loc[len(fit_df)] = fitPrntList(x, 2, vel, outputfile, model2,
                                          fitData_flip[yDataColList[4]].max(),
                                          fitData_flip['z-position(nm)'].max())
    fit_df.loc[len(fit_df)] = fitPrntList(x, 3, vel, outputfile, model3,
                                          fitData_flip[yDataColList[4]].max(),
                                          fitData_flip['z-position(nm)'].max())
    fit_df.loc[len(fit_df)] = fitPrntList(x, 'improved 1', vel, outputfile, modelB1,
                                          fitData_flip[yDataColList[4]].max(),
                                          fitData_flip['z-position(nm)'].max())
    fit_df.loc[len(fit_df)] = fitPrntList(x, 'improved 2', vel, outputfile, modelB2,
                                          fitData_flip[yDataColList[4]].max(),
                                          fitData_flip['z-position(nm)'].max())
    fit_df.loc[len(fit_df)] = fitPrntList(x, 'improved 3', vel, outputfile, modelB3,
                                          fitData_flip[yDataColList[4]].max(),
                                          fitData_flip['z-position(nm)'].max())
    fit_df.to_pickle(path.join(csvDir, "dummy.pkl"))

    plotFits(currentfile, outputfile, fitData, fitData_flip['z-position(nm)'],
             fitData_flipXY['z-position(nm)'], fitData_flip[yDataColList[4]],
             fitData_flipXY[yDataColList[4]], model1.best_fit, model2.best_fit,
             model3.best_fit, modelB1.best_fit, modelB2.best_fit, modelB3.best_fit)
    plt.savefig(path.join(imgDir, currentpic))
    # plt.show()
    plt.close()

    print(currentfile + '  Done')
