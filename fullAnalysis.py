import time
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Tkinter, tkFileDialog
from os import *
from scipy import stats
from scipy.optimize import curve_fit

def WLCmodel(x, x_off, P, L0):
    kbT = 300.0*1.380*10**(-23)
    A = 1/(4*(1-(x-x_off)/L0)**2)
    B = (x - x_off)/L0
    return kbT*(A - 0.25 + B)/P

def outputFiles(dataFiles, addon):
    L = []
    for x in range(len(dataFiles)):
        temp = dataFiles[x]
        L.append(temp[:-4] + addon)
    return L

def smooth(x,window_len,window):

    if x.ndim != 1:
        raise ValueError('smooth only accepts 1-D array')
    if x.size < window_len:
        raise ValueError('input vectoir mut be larger than window size')
    if window_len < 3:
        return quantity
    if window_len%2 == 0:
        raise ValueError('Window length must be an odd number.')
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window must be: 'flat', 'hanning', 'hamming',\
                         'bartlett', or 'blackman'")

    end = window_len - int(window_len)/2

    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]

    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')

    return y[window_len-end:-window_len+end]

## Designate input and output directories.

root = Tkinter.Tk()
root.withdraw()

info = 'Please select the folder that contains \
the data files you wish to analyze.'

srcDir = tkFileDialog.askdirectory(parent=root, initialdir="/", title=info)
dstDir = path.join(srcDir, 'images')
csvDir = path.join(srcDir, 'RetractCSVs')

## Get file list from source directory

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

##Create output files' names
dataImg = outputFiles(dataFiles, '.png')
csvOutput = outputFiles(dataFiles, '.csv')

##Pandas DataFrame to store measurements
df = pd.DataFrame(columns=['file','rupture force','location', 'WLC-P', 'WLC-L0'])
dfrow = 0

for x1 in range(len(dataFiles)):
    currentfile = dataFiles[x1]
    currentpic  = dataImg[x1]
    outputfile  = csvOutput[x1]

    dataFile = np.genfromtxt(path.join(srcDir,currentfile), skip_header=1)

    timeCh1 = dataFile[:,0]                 # s
    distance = dataFile[:,1]                # V
    timeDefl = dataFile[:,2]                # s
    deflection = dataFile[:,3]*1000000000   # nm

    ## Stiffness

    k_L = 0.07896 #N/m

    ## Find boundaries of setup, approach, retract regions
    ## using z-piezo position

    numPoints = len(distance)

    points = 500

    x_points = np.zeros(points)
    y_points = np.zeros(points)
    int_points = np.zeros(points)
    bound_pts = np.zeros((5,2))
    Vbounds = np.zeros((5,3))

    for x in range(points):
        low = x*numPoints/points
        high = low + numPoints/points
        if high > numPoints:
            high = numPoints
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            timeCh1[low:high], distance[low:high])
        x_points[x] = np.mean(timeCh1[low:high])
        y_points[x] = slope
        int_points[x] = int(slope)

    slopeSwitch = True
    bndCnt = 0
    for x in range(len(int_points)):
        if slopeSwitch:
            if int_points[x] == 0:
                bound_pts[bndCnt,0] = x_points[x]
                bound_pts[bndCnt,1] = int_points[x]
                slopeSwitch = False
                bndCnt=bndCnt+1
        else:
            if abs(int_points[x])>0:
                bound_pts[bndCnt,0] = x_points[x-1]
                bound_pts[bndCnt,1] = int_points[x-1]
                slopeSwitch = True
                bndCnt=bndCnt+1
        if bndCnt > 4:
            break
	
    if 0 in bound_pts[:,0]:
        for x in range(len(y_points)):
            int_points[x] = round(y_points[x],1)
        slopeSwitch = True
        bndCnt = 0
        for x in range(len(int_points)):
            if slopeSwitch:
                if abs(int_points[x]) < 0.2:
                    bound_pts[bndCnt,0] = x_points[x]
                    bound_pts[bndCnt,1] = int_points[x]
                    slopeSwitch = False
                    bndCnt = bndCnt + 1
            else:
                if abs(int_points[x]) > 0.2:
                    bound_pts[bndCnt,0] = x_points[x-1]
                    bound_pts[bndCnt,1] = int_points[x-1]
                    slopeSwitch = True
                    bndCnt = bndCnt + 1
            if bndCnt > 4:
                break

    bndCnt = 0
    for x in range(len(timeCh1)):
        if timeCh1[x] > bound_pts[bndCnt,0]:
            Vbounds[bndCnt,0] = timeCh1[x-1]
            Vbounds[bndCnt,1] = distance[x-1]
            Vbounds[bndCnt,2] = x
            #bndSwitch = False
            bndCnt=bndCnt+1
        if bndCnt > 4:
            break

    ## Rescaled vectors to simplify following functions
    approachT = timeCh1[Vbounds[1,2]:Vbounds[2,2]]
    approachZ = distance[Vbounds[1,2]:Vbounds[2,2]]
    approachD = deflection[Vbounds[1,2]:Vbounds[2,2]]
    retractT = timeCh1[Vbounds[3,2]:Vbounds[4,2]]
    retractZ = distance[Vbounds[3,2]:Vbounds[4,2]]
    retractD = deflection[Vbounds[3,2]:Vbounds[4,2]]

    ## Derivative of Approach curve, with rescaled Z-position
    numPoints1 = len(approachZ)
    points1 = 200
    x1_points = np.zeros(points1)
    y1_points = np.zeros(points1)
    int_points1 = np.zeros(points1)
    Abound_pts = np.zeros((4,3))
    for x in range(points1):
        low = x*numPoints1/points1
        high = low + numPoints1/points1
        if high > numPoints1:
            high = numPoints1
        s1, i1, r1, p1, e1 = stats.linregress(
            4.77*13*approachZ[low:high], approachD[low:high])
        x1_points[x] = 4.77*13*np.mean(approachZ[low:high])
        y1_points[x] = s1
        int_points1[x] = round(s1)
    
    # Find contact point, +scanning
    Abound_pts[3,0] = len(int_points1)-1
    for x in range(len(int_points1)):
        pos = x+4
        if pos > len(int_points1)-2:
            break
        if int_points1[pos] < int_points1[pos+1]:
            Abound_pts[1,0] = pos
            break
    
    # Find contact point, -scanning
    if int_points1[len(int_points1)-4] == 0:
        switch = True
    else:
        switch = False

    for x in range(len(int_points1)-5,-1,-1):
        if switch:
            if int_points1[x] < int_points1[x-1]:
                Abound_pts[3,0] = x-1
                switch = False
        else:
            if int_points1[x] < 1:
                Abound_pts[2,0] = x+1
                break

    for x in range(len(Abound_pts)):
        Abound_pts[x,1] = x1_points[Abound_pts[x,0]]
        Abound_pts[x,2] = int_points1[Abound_pts[x,0]]

    ## Derivative of Retract curve, with rescaled Z-position
    numPoints2 = len(retractZ)
    points2 = 200
    x2_points = np.zeros(points2)
    y2_points = np.zeros(points2)
    int_points2 = np.zeros(points2)
    Rbound_pts = np.zeros((4,3))
    for x in range(points2):
        low = x*numPoints2/points2
        high = low + numPoints2/points2
        if high > numPoints2:
            high = numPoints2
        s2, i2, r2, p2, e2 = stats.linregress(4.77*13*retractZ[low:high], retractD[low:high])
        x2_points[x] = 4.77*13*np.mean(retractZ[low:high])
        y2_points[x] = s2
        int_points2[x] = round(s2)
    
    # Find contact point, +scanning
    Rbound_pts[0,0] = len(int_points2)-1
    Rbound_pts[3,0] = 0
    if int(np.average(int_points2[:4])) < 1:
        switch = True
    else:
        switch = False

    for x in range(len(int_points2)):
        pos = x+4
        if pos > len(int_points2)-3:
            break
        if switch:
            if (int_points2[pos] < int_points2[pos+1]) and (int_points2[pos] < int_points2[pos+2]):
                Rbound_pts[3,0] = pos+1
                switch = False
        else:
            if int_points2[pos] < 1:
                Rbound_pts[2,0] = pos-1
                break

    # Find contact point, -scanning
    for x in range(len(int_points2)-4,-1,-1):
        if (int_points2[x] < int_points2[x-1]): #and (int_points2[x-1] > 0)
            Rbound_pts[1,0] = x
            break

    for x in range(len(Rbound_pts)):
        Rbound_pts[x,1] = x2_points[Rbound_pts[x,0]]
        Rbound_pts[x,2] = int_points2[Rbound_pts[x,0]]

    Rbounds = np.zeros((4,3))
    rBndCnt = 0
    for x in range(len(retractZ)-1,-1,-1):
        if 4.77*13*retractZ[x-1] > Rbound_pts[rBndCnt,1]:
            Rbounds[rBndCnt,0] = x
            Rbounds[rBndCnt,1] = 4.77*13*retractZ[x]
            Rbounds[rBndCnt,2] = retractD[x]
            rBndCnt = rBndCnt + 1
        if rBndCnt > 3:
            break

    # Find baseline & contact region linear regression
    baselineS, baselineI, baselineR, baselineP, baselineE = stats.linregress(
            4.77*13*retractZ[Rbounds[1,0]:Rbounds[0,0]], retractD[Rbounds[1,0]:Rbounds[0,0]])
    contactS, contactI, contactR, contactP, contactE = stats.linregress(
            4.77*13*retractZ[Rbounds[3,0]:Rbounds[2,0]], retractD[Rbounds[3,0]:Rbounds[2,0]])

    x_shift = (baselineI - contactI)/(contactS - baselineS)
    y_shift = contactS * x_shift + contactI

    ## Linear Regression on approach/retract regions
    ## __1 = slope ; __2 = intercept ; __3 = r_value ; __4 = p_value ; __5 = std_error
    ##
    ## Setup Speed
    setup1,setup2,setup3,setup4,setup5 = stats.linregress(
        timeCh1[0:Vbounds[0,2]], distance[0:Vbounds[0,2]])
    #print "setup v =", abs(setup1*13*4.77), "nm/s"

    ## Approach Speed
    appr1,appr2,appr3,appr4,appr5 = stats.linregress(
        approachT, approachZ)
    #print "approch v =", abs(appr1*13*4.77), "nm/s"

    ## Retract Speed
    retr1,retr2,retr3,retr4,retr5 = stats.linregress(
        retractT, retractZ)
    #print "retract v =", abs(retr1*13*4.77), "nm/s"
    
    smooth11 = smooth((k_L*(retractD - y_shift)/contactS),11,'hanning')
    smooth25 = smooth((k_L*(retractD - y_shift)/contactS),25,'hanning')
    smooth55 = smooth((k_L*(retractD - y_shift)/contactS),55,'hanning')
    smooth75 = smooth((k_L*(retractD - y_shift)/contactS),75,'hanning')
    
    # Find Rupture Force	
    ruptureI = np.argmin(retractD[Rbounds[3,0]:])
    ruptureF = k_L*(retractD[ruptureI] - y_shift)/contactS
    ruptureL = (4.77*13*retractZ[ruptureI] - (retractD[ruptureI] - y_shift) - x_shift)
    
    # Fit WLC model to rupture
##    fitXdata = (4.77*13*retractZ[(ruptureI-200):ruptureI] - (retractD[ruptureI] - y_shift) - x_shift)
##    fitYdata = k_L*(retractD[(ruptureI-200):ruptureI] - y_shift)/contactS
##    popt, pcov = curve_fit(WLCmodel, fitXdata, fitYdata)
##
##    print 'popt =', popt
    
    WLC_P = 1.0
    WLC_L0 = 2.0
    
    #Add data to pandas DataFrame
    df.loc[dfrow] = [currentfile, ruptureF, ruptureL, WLC_P, WLC_L0]
    dfrow = dfrow + 1

    #Output Calculations
    output = np.column_stack((4.77*13*retractZ - x_shift,
                             (4.77*13*retractZ - (retractD - y_shift) - x_shift),
                              k_L*(retractD - y_shift)/contactS,
                              smooth11,smooth25,smooth55,smooth75))

    csvheader = "z-position(nm),separation(nm),Force(nN),Force_11(nN),Force_25(nN),Force_55(nN),Force_75(nN),v=%d nm/s" % (abs(retr1*13*4.77))

    np.savetxt(path.join(csvDir,outputfile), output, header=csvheader, comments="", delimiter=',')

    ## Figures

    plt.figure(figsize=(20,10))
    plt.subplot(2,3,2)
    plt.title("Approach Slope")
    plt.plot(x1_points, y1_points, 'b.')
    plt.plot(x1_points, int_points1, 'r.')
    plt.plot(Abound_pts[:,1],Abound_pts[:,2], 'go')
    plt.ylabel("Defl. Deriv")
    plt.xlabel("Z-Position (nm)")

    plt.subplot(2,3,3)
    plt.title("Retract Slope")
    plt.plot(x2_points, y2_points, 'b.')
    plt.plot(x2_points, int_points2, 'r.')
    plt.plot(Rbound_pts[:,1],Rbound_pts[:,2], 'go')
    plt.ylabel("Defl. Deriv")
    plt.xlabel("Z-Position (nm)")

    plt.subplot(2,3,1)
    plt.title("Z-position (V)")
    plt.plot(timeCh1, distance)
    plt.plot(Vbounds[:,0], Vbounds[:,1], 'ro')
    plt.xlabel("Time (s)")
    plt.ylabel("Z-Position (V)")

    plt.subplot(2,3,6)
    plt.title("Retract")
    plt.plot(4.77*13*retractZ,retractD)
    plt.plot(4.77*13*retractZ - x_shift,retractD - y_shift)
    plt.plot(0, 0, 'ro')
    plt.plot(Rbounds[:,1], Rbounds[:,2], 'ro')
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (V)")

    plt.subplot(2,3,5)
    plt.title("Approach")
    plt.plot(approachZ,approachD)
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (V)")

    plt.subplot(2,3,4)
    plt.title("Setup")
    plt.plot(distance[0:Vbounds[0,2]], deflection[0:Vbounds[0,2]])
    plt.ylabel("Deflection (nm)")
    plt.xlabel("Z-position (V)")

    plt.savefig(path.join(dstDir,currentpic))
    plt.close()
	
    print "Completed ", x1+1, " of ", len(dataFiles), " files."

df.to_excel(path.join(csvDir,'dataframe.xlsx'), sheet_name='Sheet1')
print "Finished analyzing", path.split(srcDir)[1]
print 'It took {:.2f} seconds to analyze %d files.'.format(time.time()-start) % (len(dataFiles))
