import sys
import time
import Tkinter, tkFileDialog
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import multiLinReg as MLR
from os import *
from scipy import stats
from scipy.optimize import curve_fit
from lmfit import Model
from polymerModels import WLCmodel, FJCmodel

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
csvRupture = outputFiles(dataFiles, '-rupture.csv')

##Pandas DataFrame to store measurements
df = pd.DataFrame(columns=['file','rupture force (pN)','location', 'speed (nm/s)', 'WLC-P', 'WLC-L0', 'x_off'])
dfrow = 0

for x1 in range(len(dataFiles)):
	currentfile = dataFiles[x1]
	currentpic  = dataImg[x1]
	outputfile  = csvOutput[x1]
	ruptureFile = csvRupture[x1]

	dataFile = np.genfromtxt(path.join(srcDir,currentfile), skip_header=1)

	timeCh1 = dataFile[:,0]                 # s
	distance = dataFile[:,1]                # V
	timeDefl = dataFile[:,2]                # s
	deflection = dataFile[:,3]*1000000000   # nm

	## Stiffness

	k_L = 0.034 #N/m

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
	approachZ = 4.77*13.0*distance[Vbounds[1,2]:Vbounds[2,2]]
	approachD = deflection[Vbounds[1,2]:Vbounds[2,2]]
	approachF = k_L * approachD
	retractT = timeCh1[Vbounds[3,2]:Vbounds[4,2]]
	retractZ = 4.77*13.0*distance[Vbounds[3,2]:Vbounds[4,2]]
	retractD = deflection[Vbounds[3,2]:Vbounds[4,2]]
	retractF = k_L * retractD

	try:
		contactS, contactI, baselineS, baselineI = MLR.multiLinReg(retractZ,retractD)
	except Exception:
		sys.exc_clear()
		print "File %s failed" % (currentfile)
		plt.plot(retractZ,retractD)
		plt.xlabel("Z-position (nm)")
		plt.ylabel("Deflection (nm)")
		plt.savefig(path.join(dstDir,currentpic))
		plt.close()
		continue

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
	gmod.set_param_hint('L_C', value = -60.0)
	gmod.set_param_hint('L_P', value = -0.38, min=-0.42, max=-0.34)
	gmod.set_param_hint('a', value = 0.0, min=-10.0, max=10.0)
	gmod.set_param_hint('b', value = 0.0, min=-10.0, max=10.0)
	params = gmod.make_params()
	try:
		result = gmod.fit(smooth25[originPt:ruptureI], x=separation[originPt:ruptureI]) # method='cobyla'
	except Exception:
		skipPLT5 = False
		sys.exc_clear()
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
	# FJCmod.set_param_hint('L0') #, value = -56.0)
	# FJCmod.set_param_hint('b') #, value = -3.8, min=-4.0, max=-3.6)
	# ##FJCmod.set_param_hint('a', value = 0.0, min=-5.0, max=5.0)
	# FJCparams = FJCmod.make_params()
	# try:
		# FJCresult = FJCmod.fit(separation[originPt:ruptureI], x=smooth25[originPt:ruptureI]) # method='cobyla'
	# except Exception:
		# print "FJC failed"
		# skipPLT6 = False
		# sys.exc_clear()
	# if skipPLT6:
		# #x_off = result.params['a'].value
		# FJC_L0 = FJCresult.params['L0'].value
		# FJC_b = FJCresult.params['b'].value
	# else:
		# #x_off = 0.0
		# FJC_L0 = 0.0
		# FJC_b = 0.0

	#Add data to pandas DataFrame
	df.loc[dfrow] = [currentfile, 1000.0*abs(ruptureF), ruptureL, abs(retr1), WLC_P, WLC_L0, x_off]
	dfrow = dfrow + 1

	#Output Calculations
	output = np.column_stack((retractZ - x_shift, separation,
							 k_L*(retractD - y_shift)/contactS,
							 smooth11,smooth25,smooth55,smooth75))
	ruptureOut = np.column_stack((separation[originPt:ruptureI],smooth25[originPt:ruptureI]))

	csvheader = "z-position(nm),separation(nm),Force(nN),Force_11(nN),Force_25(nN),Force_55(nN),Force_75(nN),v=%d nm/s" % (abs(retr1))
	
	ruptureH = "separation(nm),Force_25(nN)"

	np.savetxt(path.join(csvDir,outputfile), output, header=csvheader, comments="", delimiter=',')
	np.savetxt(path.join(csvDir,ruptureFile), ruptureOut, header=ruptureH, comments="", delimiter=',')

	## Figures

	plt.figure(figsize=(20,10))
	plt.subplot(2,3,1)
	plt.title("Z-position (nm)")
	plt.plot(timeCh1, distance)
	plt.plot(Vbounds[:,0], Vbounds[:,1], 'ro')
	plt.xlabel("Time (s)")
	plt.ylabel("Z-Position (V)")
	
	plt.subplot(2,3,2)
	plt.title("Fitted Retract")
	plt.plot(retractZ,retractD)
	plt.plot(retractZ,baselineS*retractZ + baselineI)
	plt.plot(retractZ,contactS*retractZ + contactI)
	plt.ylabel("Deflection (nm)")
	plt.xlabel("Z-position (nm)")
	plt.axis([min(retractZ)-5,max(retractZ)+5,min(retractD)-10,max(retractD)+10])
	
	plt.subplot(2,3,3)
	plt.title("Full Retract")
	plt.plot(retractZ,retractD)
	plt.plot(retractZ - x_shift,retractD - y_shift)
	plt.plot(0, 0, 'ro')
	plt.ylabel("Deflection (nm)")
	plt.xlabel("Z-position (nm)") 
	
	plt.subplot(2,3,4)
	plt.title("Retract")
	plt.plot(retractZ - x_shift,retractD - y_shift)
	plt.plot(0, 0, 'ro')
	plt.ylabel("Deflection (nm)")
	plt.xlabel("Z-position (nm)") 
	plt.axis([-100,10,min(retractD - y_shift)-5,30])

	if skipPLT5:
		plt.subplot(2,3,5)
		plt.title("Fit")
		plt.plot(separation[originPt:ruptureI], smooth25[originPt:ruptureI], 'b.')
		plt.plot(separation[originPt:ruptureI], result.init_fit, 'k--')
		plt.plot(separation[originPt:ruptureI], result.best_fit, 'r-')
		plt.ylabel("Force (nN)")
		plt.xlabel("Separation (nm)")
	else:
		plt.subplot(2,3,5)
		plt.title("Fit")
		plt.plot(separation[originPt:ruptureI], smooth25[originPt:ruptureI], 'b.')
		plt.ylabel("Force (nN)")
		plt.xlabel("Separation (nm)")
	
	# if skipPLT6:
		# plt.subplot(2,3,6)
		# plt.title("Fit")
		# plt.plot(separation[originPt:ruptureI], smooth25[originPt:ruptureI], 'b.')
		# #plt.plot(separation[originPt:ruptureI], FJCresult.init_fit, 'k--')
		# plt.plot(FJCresult.best_fit, smooth25[originPt:ruptureI], 'r-')
		# plt.ylabel("Force (nN)")
		# plt.xlabel("Separation (nm)")
	# else:
		# plt.subplot(2,3,6)
		# plt.title("Fit")
		# plt.plot(separation[originPt:ruptureI], smooth25[originPt:ruptureI], 'b.')
		# plt.ylabel("Force (nN)")
		# plt.xlabel("Separation (nm)")
	
	plt.savefig(path.join(dstDir,currentpic))
	plt.close()

	print "Completed ", x1+1, " of ", len(dataFiles), " files."

df.to_excel(path.join(csvDir,'dataframe.xlsx'), sheet_name='Sheet1')
print "Finished analyzing", path.split(srcDir)[1]
print 'It took {:.2f} seconds to analyze %d files.'.format(time.time()-start) % (len(dataFiles))
