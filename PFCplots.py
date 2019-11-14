import matplotlib.pyplot as plt


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
             model1fit, model2fit, model3fit, modelB1fit, modelB2fit, modelB3fit):
    """Info
    """
    plt.figure(figsize=(12, 16)).suptitle(outputfile[:-4] + '  '
                                         + fitData.columns[-1])
    plt.title(currentfile)
    plt.subplot(3, 2, 1)
    plt.plot(xdata, ydata, 'b,')
    plt.plot(xdata, model1fit, 'r-', label='best fit')
    plt.axis([0, max(xdata)+5, 0, max(ydata)*1.1])
    plt.legend(loc='best')

    plt.subplot(3, 2, 3)
    plt.plot(xdata, ydataXY, 'b,')
    plt.plot(xdata, model2fit, 'r-', label='best fit')
    plt.axis([0, max(xdata)+5, 0, max(ydataXY)*1.1])
    plt.legend(loc='best')

    plt.subplot(3, 2, 5)
    plt.plot(xdataXY, ydataXY, 'b,')
    plt.plot(xdataXY, model3fit, 'r-', label='best fit')
    plt.axis([0, max(xdataXY)+5, 0, max(ydataXY)*1.1])
    plt.legend(loc='best')

    plt.subplot(3, 2, 2)
    plt.plot(xdata, ydata, 'b,')
    plt.plot(xdata, modelB1fit, 'r-', label='best fit')
    plt.axis([0, max(xdata)+5, 0, max(ydata)*1.1])
    plt.legend(loc='best')

    plt.subplot(3, 2, 4)
    plt.plot(xdata, ydataXY, 'b,')
    plt.plot(xdata, modelB2fit, 'r-', label='best fit')
    plt.axis([0, max(xdata)+5, 0, max(ydataXY)*1.1])
    plt.legend(loc='best')

    plt.subplot(3, 2, 6)
    plt.plot(xdataXY, ydataXY, 'b,')
    plt.plot(xdataXY, modelB3fit, 'r-', label='best fit')
    plt.axis([0, max(xdataXY)+5, 0, max(ydataXY)*1.1])
    plt.legend(loc='best')
