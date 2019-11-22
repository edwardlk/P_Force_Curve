from multiLinReg import multiLinRegDirect
from os import path
import pandas as pd
import numpy as np
from PFCfuncs import smooth
from PFCplots import plotEverything
import matplotlib.pyplot as plt

fileTest = True

k_L = 0.034
print('k_L = {}'.format(k_L))

if fileTest:
    srcDir = R'D:\mlr_fix_test'

imgDir = path.join(srcDir, 'images')
csvDir = path.join(srcDir, 'RetractCSVs')

outputCSV = path.join(csvDir, 'dataframe.xlsx')
outputPkl = path.join(csvDir, "dummy.pkl")

outputDF = pd.read_excel(outputCSV, index_col=0)
outputDF.to_pickle(outputPkl)

for x in range(len(outputDF)):
    testCB = (outputDF['cStart'][x] == outputDF['bStart'][x])
    if not testCB:
        dataFile = np.genfromtxt(path.join(srcDir, outputDF['filename'][x]),
                                 skip_header=1)

        (timeCh1, distance, timeDefl, deflection) = dataFile.T
        deflection = deflection*1000000000  # convert deflection to nm

        # outputDF['Vb3'][x]:outputDF['Vb4'][x]
        VboundsI = [outputDF['Vb1'][x], outputDF['Vb2'][x], outputDF['Vb3'][x],
                    outputDF['Vb4'][x], outputDF['Vb5'][x]]
        print(VboundsI)
        VboundsXY = np.zeros((5, 2))
        for x1 in range(len(VboundsI)):
            VboundsXY[x1, 0] = timeCh1[VboundsI[x1]]
            VboundsXY[x1, 1] = distance[VboundsI[x1]]

        retractT = timeCh1[outputDF['Vb4'][x]:outputDF['Vb5'][x]]
        retractZ_orig = 4.77*13*distance[outputDF['Vb4'][x]:outputDF['Vb5'][x]]
        retractD_orig = deflection[outputDF['Vb4'][x]:outputDF['Vb5'][x]]

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

        (contactS, contactI, baselineS, baselineI) = multiLinRegDirect(
            retractZ_orig, retractD_orig,
            outputDF['bStart'][x], outputDF['cStart'][x])

        x_shift = (baselineI - contactI)/(contactS - baselineS)
        y_shift = contactS * x_shift + contactI

        retractZ = retractZ_orig - x_shift
        retractD = retractD_orig - y_shift

        originPt = abs(retractZ).argmin()
        retractD = (retractD - baselineS * retractZ) / (contactS - baselineS)

        separation = retractZ - retractD

        # Retract Speed
        retr1 = 0.0
        # retr1, retr2, retr3, retr4, retr5 = stats.linregress(retractT, retractZ)

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

        # Figures
        plotEverything('currentpic', abs(retr1), originPt, baselineS, baselineI,
                       contactS, contactI, retractZ_orig, retractD_orig, smooth3,
                       separation, timeCh1, distance, retractZ, retractD,
                       VboundsXY, VboundsI)
        plt.show()
        # plt.savefig(path.join(dstDir, currentpic))
        # plt.close()
