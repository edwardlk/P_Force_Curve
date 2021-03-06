import matplotlib.pyplot as plt
import numpy as np


def plotEverything(currentpic, v, originPt, baselineS, baselineI,
                   contactS, contactI, retractZ_orig, retractD_orig, smooth3,
                   separation, timeCh1, distance, retractZ, retractD,
                   VboundsXY, VboundsI):
    """ Info
    """

    plt.figure(figsize=(20, 10))
    plt.suptitle('{} - v = {:.1f} m/s'.format(currentpic[:-11], v),
                 fontsize=24)
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
    try:
        plt.axis([-150, 10, min(retractD)-5, 20])
    except ValueError:
        plt.axis([-100, 10, -5, 20])
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(15))
    plt.grid(True, which="both")

    plt.subplot(2, 3, 5)
    plt.title("Retract")
    plt.plot(retractZ, retractD, ',')
    plt.plot(0, 0, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Separation (nm)")
    try:
        plt.axis([-100, 10, min(retractD)-5, 20])
    except ValueError:
        plt.axis([-100, 10, -5, 20])
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(10))
    plt.grid(True, which="both")

    minS3 = min(smooth3)
    maxS3 = max(smooth3)
    plt.subplot(2, 3, 6)
    plt.title("Retract")
    # plt.plot(separation, smooth3 + abs(minS3*1.1), ',b')
    plt.plot(retractZ, smooth3, ',k')
    plt.plot(0, 0, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Smoothed Force (nN)")
    try:
        plt.axis([-50, 10, minS3*1.1, maxS3*0.1])
    except ValueError:
        plt.axis([-50, 10, -5, 10])
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(4))
    plt.grid(True, which="both")


def plotEverything2(currentpic, v, originPt, baselineS, baselineI,
                    contactS, contactI, retractZ_orig, retractD_orig, smooth3,
                    separation, distance, retractZ, retractD, x_shift, y_shift):
    """
    """
    plt.figure(figsize=(20, 10))
    plt.suptitle('{}'.format(currentpic[:-4]), fontsize=24)
    grid = plt.GridSpec(2, 3)
    plt.subplot(grid[0, 0])
    plt.title("Fitted Retract")
    plt.plot(retractZ_orig, retractD_orig, ',')
    plt.plot(retractZ_orig, baselineS*retractZ_orig + baselineI)
    plt.plot(retractZ_orig, contactS*retractZ_orig + contactI)
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (nm)")
    plt.axis([min(retractZ_orig)-5, max(retractZ_orig)+5,
              min(retractD_orig)-10, max(retractD_orig)+10])
    plt.grid(True, which="both")

    plt.subplot(grid[0, 1])
    plt.title("Fitted Retract")
    plt.plot(retractZ_orig, retractD_orig, ',')
    plt.plot(retractZ_orig, baselineS*retractZ_orig + baselineI)
    plt.plot(retractZ_orig, contactS*retractZ_orig + contactI)
    plt.plot(x_shift, y_shift, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (nm)")
    plt.axis([x_shift-15, x_shift+150, y_shift-10, y_shift+10])
    plt.grid(True, which="both")

    plt.subplot(grid[0, 2])
    plt.title("Shifted Retract")
    plt.plot(retractZ_orig, retractD_orig, ',')
    plt.plot(retractZ, retractD, ',')
    plt.plot(0, 0, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (nm)")
    plt.grid(True, which="both")

    temp = abs(retractZ - 16.0)
    maxID = np.argmin(temp)

    plt.subplot(grid[1, 0:])
    plt.title("Retract")
    plt.plot(retractZ, retractD, ',')
    plt.plot(0, 0, 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (nm)")
    try:
        plt.axis([-160, 16, min(retractD)-3, retractD[maxID]])
    except ValueError:
        plt.axis([-160, 16, -5, 15])
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(4))
    plt.grid(True, which="both")

    # plt.subplot(2, 3, 5)
    # plt.title("Retract")
    # plt.plot(retractZ, retractD, ',')
    # plt.plot(0, 0, 'ro')
    # plt.ylabel("Deflection (nm)")
    # plt.xlabel("Z-position (nm)")
    # try:
    #     plt.axis([-100, 10, min(retractD)-5, 20])
    # except ValueError:
    #     plt.axis([-100, 10, -5, 20])
    # plt.gca().xaxis.set_major_locator(plt.MultipleLocator(10))
    # plt.grid(True, which="both")
    #
    # minS3 = min(smooth3)
    # maxS3 = max(smooth3)
    # plt.subplot(2, 3, 6)
    # plt.title("Retract")
    # # plt.plot(separation, smooth3 + abs(minS3*1.1), ',b')
    # plt.plot(retractZ, smooth3, ',k')
    # plt.plot(0, 0, 'ro')
    # plt.ylabel("Smoothed Force (nN)")
    # plt.xlabel("Z-position (nm)")
    # try:
    #     plt.axis([-52, 12, minS3*1.5, maxS3*0.1])
    # except ValueError:
    #     plt.axis([-52, 12, -5, 10])
    # plt.gca().xaxis.set_major_locator(plt.MultipleLocator(4))
    # plt.grid(True, which="both")


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
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(25))
    plt.grid(True, which="both")


def plotFits(currentfile, outputfile, fitData, xdata, xdataXY, ydata, ydataXY,
             model1fit, model2fit, model3fit, modelB1fit, modelB2fit,
             modelB3fit, x_data, y_data,
             minGuessID, minGuessRange, fitStartID, minID):
    """Info
    """
    vel = float(fitData.columns[-2][2:-5])
    title1 = (outputfile[:-4] + '  v = {:.2f} '
              + fitData.columns[-1]).format(vel)

    plt.figure(figsize=(12, 15)).suptitle(title1, fontsize=24)
    plt.title(currentfile)
    plt.subplot(411)
    plt.plot(x_data, y_data, 'b,')
    plt.plot(x_data.iloc[[minGuessID]], min(y_data)*1.1, 'rx')
    plt.plot(x_data.iloc[[minGuessID+minGuessRange]], min(y_data)*1.1, 'r|')
    plt.plot(x_data.iloc[[minGuessID-minGuessRange]], min(y_data)*1.1, 'r|')
    plt.plot(x_data.iloc[[minID]], min(y_data)*1.1, 'gx')
    plt.plot(x_data.iloc[[fitStartID]], min(y_data)*1.1, 'g|')
    plt.plot(x_data[fitStartID:minID], y_data[fitStartID:minID], 'r,')
    plt.axis([(x_data.iloc[[minID]][minID]-10), 5,
              min(y_data)*1.2, abs(min(y_data))])
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(2))
    y_maj = abs(y_data.at[minID])
    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(y_maj))
    plt.grid(True, axis='both')

    plt.subplot(423)
    plt.plot(xdata, ydata, 'b,')
    plt.plot(xdata, model1fit, 'r-', label='best fit')
    plt.axis([0, max(xdata)+5, 0, max(ydata)*1.1])
    plt.legend(loc='best')

    plt.subplot(425)
    plt.plot(xdata, ydataXY, 'b,')
    plt.plot(xdata, model2fit, 'r-', label='best fit')
    plt.axis([0, max(xdata)+5, 0, max(ydataXY)*1.1])
    plt.legend(loc='best')

    plt.subplot(427)
    plt.plot(xdataXY, ydataXY, 'b,')
    plt.plot(xdataXY, model3fit, 'r-', label='best fit')
    plt.axis([0, max(xdataXY)+5, 0, max(ydataXY)*1.1])
    plt.legend(loc='best')

    plt.subplot(424)
    plt.plot(xdata, ydata, 'b,')
    plt.plot(xdata, modelB1fit, 'r-', label='best fit')
    plt.axis([0, max(xdata)+5, 0, max(ydata)*1.1])
    plt.legend(loc='best')

    plt.subplot(426)
    plt.plot(xdata, ydataXY, 'b,')
    plt.plot(xdata, modelB2fit, 'r-', label='best fit')
    plt.axis([0, max(xdata)+5, 0, max(ydataXY)*1.1])
    plt.legend(loc='best')

    plt.subplot(428)
    plt.plot(xdataXY, ydataXY, 'b,')
    plt.plot(xdataXY, modelB3fit, 'r-', label='best fit')
    plt.axis([0, max(xdataXY)+5, 0, max(ydataXY)*1.1])
    plt.legend(loc='best')
