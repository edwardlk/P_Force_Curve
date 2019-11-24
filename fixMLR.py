from multiLinReg import multiLinRegDirect
from os import path
import pandas as pd
import numpy as np
from PFCfuncs import smooth
from PFCplots import plotEverything
import matplotlib.pyplot as plt

fileTest = True

# k_L = 0.034
# print('k_L = {}'.format(k_L))

if fileTest:
    srcDir = R'F:\_data\Avidin-Biotin\fixmlrtest'

imgDir = path.join(srcDir, 'images')
csvDir = path.join(srcDir, 'RetractCSVs')

outputCSV = path.join(csvDir, 'dataframe.xlsx')
outputPkl = path.join(csvDir, "dummy.pkl")

outputDF = pd.read_excel(outputCSV, index_col=0)
outputDF.to_pickle(outputPkl)

# 0     1           2               3           4       5
# f_num filename    min_location    fit_start   cStart	bStart
# 6     7   8   9   10
# Vb1   Vb2 Vb3 Vb4 Vb5

for x in range(len(outputDF)):
    f_name = outputDF.iloc[x, 1]
    min_loc = outputDF.iloc[x, 2]
    fit_st = outputDF.iloc[x, 3]
    cstart = outputDF.iloc[x, 4]
    bstart = outputDF.iloc[x, 5]
    VboundsI = [outputDF.iloc[x, 6], outputDF.iloc[x, 7], outputDF.iloc[x, 8],
                outputDF.iloc[x, 9], outputDF.iloc[x, 10]]

    if not cstart == bstart:
        print('Fixing {}...'.format(f_name), end='')

        currentpic = f_name[:-4] + '.png'
        outputData = f_name[:-4] + '.csv'

        dataFile = np.genfromtxt(path.join(srcDir, f_name), skip_header=1)

        # [8]: v=31.161730724995547 nm/s [9]: k_L=0.034
        oldCSV = pd.read_csv(path.join(csvDir, outputData))
        oldCSVh = oldCSV.columns

        retr1 = float(oldCSVh[8][2:-5])
        k_L = float(oldCSVh[9][4:])

        (timeCh1, distance, timeDefl, deflection) = dataFile.T
        deflection = deflection*1000000000  # convert deflection to nm

        VboundsXY = np.zeros((5, 2))
        for x1 in range(len(VboundsI)):
            VboundsXY[x1, 0] = timeCh1[VboundsI[x1]]
            VboundsXY[x1, 1] = distance[VboundsI[x1]]

        retractT = timeCh1[VboundsI[3]:VboundsI[4]]
        retractZ_orig = 4.77*13*distance[VboundsI[3]:VboundsI[4]]
        retractD_orig = deflection[VboundsI[3]:VboundsI[4]]

        # Remove data if deflection out of range
        maxDefl = max(retractD_orig)
        x = 0
        while retractD_orig[x1] == maxDefl:
            x1 += 1
            if x1 > len(retractD_orig)-1:
                break

        retractT = retractT[x1:]
        retractZ_orig = retractZ_orig[x1:]
        retractD_orig = retractD_orig[x1:]

        (contactS, contactI, baselineS, baselineI) = multiLinRegDirect(
            retractZ_orig, retractD_orig, cstart, bstart)

        x_shift = (baselineI - contactI)/(contactS - baselineS)
        y_shift = contactS * x_shift + contactI

        retractZ = retractZ_orig - x_shift
        retractD = retractD_orig - y_shift

        originPt = abs(retractZ).argmin()
        retractD = (retractD - baselineS * retractZ) / (contactS - baselineS)

        separation = retractZ - retractD

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

        output = np.column_stack((retractZ, separation, retractD, k_L*retractD,
                                  smooth1, smooth2, smooth3, smooth4))
        csvheader = ('z-position(nm),separation(nm),retractD,'
                     'Force(nN),Force_{}(nN),Force_{}(nN),Force_{}(nN),'
                     'Force_{}(nN),v={} nm/s,k_L={}')
        csvheader = csvheader.format(
            smDict['smooth1'], smDict['smooth2'], smDict['smooth3'],
            smDict['smooth4'], abs(retr1), k_L)

        np.savetxt(path.join(csvDir, outputData), output, header=csvheader,
                   comments="", delimiter=',')

        # Figures
        plotEverything(currentpic, abs(retr1), originPt, baselineS, baselineI,
                       contactS, contactI, retractZ_orig, retractD_orig,
                       smooth3, separation, timeCh1, distance, retractZ,
                       retractD, VboundsXY, VboundsI)
        # plt.show()
        # print(path.join(imgDir, currentpic))
        plt.savefig(path.join(imgDir, currentpic))
        plt.close()

        print('Done')
