import sys
import time
from tkinter import Tk, filedialog
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import multiLinReg as MLR
from os import path, listdir, makedirs
from scipy import stats
from scipy.optimize import curve_fit
from lmfit import Model
from polymerModels import WLCmodel, FJCmodel
from PFCfuncs import outputFiles, smooth, returnBoundaries

# Designate input and output directories.
root = Tk()
root.withdraw()

info = ("Please select the folder that contains the data files "
        "you wish to analyze.")

srcDir = filedialog.askdirectory(parent=root, initialdir="/", title=info)
dstDir = path.join(srcDir, 'images')
csvDir = path.join(srcDir, 'RetractCSVs')

# Get file list from source directory

dataFiles = listdir(srcDir)
dataFiles.sort()

start = time.time()

if not path.exists(dstDir):
    makedirs(dstDir)
else:
    dataFiles.remove('images')
if not path.exists(csvDir):
    makedirs(csvDir)
else:
    dataFiles.remove('RetractCSVs')

# Create output files' names
dataImg = outputFiles(dataFiles, '.png')
csvOutput = outputFiles(dataFiles, '.csv')
csvRupture = outputFiles(dataFiles, '-rupture.csv')

# Pandas DataFrame to store measurements
df = pd.DataFrame(columns=['file', 'rupture force (pN)', 'location',
                           'speed (nm/s)', 'WLC-P', 'WLC-L0', 'x_off'])
dfrow = 0

for x1 in range(len(dataFiles)):
    currentfile = dataFiles[x1]
    currentpic = dataImg[x1]
    outputfile = csvOutput[x1]
    ruptureFile = csvRupture[x1]

    dataFile = np.genfromtxt(path.join(srcDir, currentfile), skip_header=1)

    (timeCh1, distance, timeDefl, deflection) = dataFile.T
    deflection = deflection*1000000000

    # Stiffness
    k_L = 0.034  # N/m

    # Find boundaries of setup, approach, retract regions
    # using z-piezo position
    Vbounds = returnBoundaries(timeCh1, distance, 500)

    # Rescaled vectors to simplify following functions
    approachT = timeCh1[int(Vbounds[1, 2]):int(Vbounds[2, 2])]
    approachZ = 4.77*13*distance[int(Vbounds[1, 2]):int(Vbounds[2, 2])]
    approachD = deflection[int(Vbounds[1, 2]):int(Vbounds[2, 2])]
    approachF = k_L * approachD
    retractT = timeCh1[int(Vbounds[3, 2]):int(Vbounds[4, 2])]
    retractZ = 4.77*13*distance[int(Vbounds[3, 2]):int(Vbounds[4, 2])]
    retractD = deflection[int(Vbounds[3, 2]):int(Vbounds[4, 2])]
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
        continue

    x_shift = (baselineI - contactI)/(contactS - baselineS)
    y_shift = contactS * x_shift + contactI

    # Linear Regression on approach/retract regions
    # __1 = slope ; __2 = intercept ; __3 = r_value ;
    # __4 = p_value ; __5 = std_error
    # Setup Speed
    setup1, setup2, setup3, setup4, setup5 = stats.linregress(
        timeCh1[0:int(Vbounds[0, 2])], distance[0:int(Vbounds[0, 2])])
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
    df.loc[dfrow] = [currentfile, 1000.0*abs(ruptureF), ruptureL, abs(retr1),
                     WLC_P, WLC_L0, x_off]
    dfrow = dfrow + 1

    # Output Calculations
    output = np.column_stack((retractZ - x_shift, separation, retractD,
                              k_L*(retractD - y_shift)/contactS, smooth11,
                              smooth25, smooth55, smooth75))
    ruptureOut = np.column_stack((separation[originPt:ruptureI],
                                  smooth25[originPt:ruptureI]))

    csvheader = ("z-position(nm),separation(nm),retractD,Force(nN),Force_11(nN),"
                 "Force_25(nN),Force_55(nN),Force_75(nN),v=%d nm/s"
                 % (abs(retr1)))

    ruptureH = "separation(nm),Force_25(nN)"

    np.savetxt(path.join(csvDir, outputfile), output, header=csvheader,
               comments="", delimiter=',')
    np.savetxt(path.join(csvDir, ruptureFile), ruptureOut, header=ruptureH,
               comments="", delimiter=',')

    # Figures

    plt.figure(figsize=(20, 10))
    plt.subplot(2, 3, 1)
    plt.title("Z-position (nm)")
    plt.plot(timeCh1, distance)
    plt.plot(Vbounds[:, 0], Vbounds[:, 1], 'ro')
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

    plt.savefig(path.join(dstDir, currentpic))
    plt.close()

    print("Completed ", x1+1, " of ", len(dataFiles), " files.")

df.to_excel(path.join(csvDir, 'dataframe.xlsx'), sheet_name='Sheet1')
print("Finished analyzing", path.split(srcDir)[1])
print('It took {:.2f} seconds to analyze %d files.'.format(time.time()-start)
      % (len(dataFiles)))
