import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import Tkinter, tkFileDialog
from os import *
from scipy import stats
from scipy.optimize import curve_fit

def multiLinReg(x_data,y_data):
	print "Running MLR"
	print len(x_data)
	midlen = len(x_data)/20000
	init_points = np.zeros((19999,4))
	for x in range(19999):
		init_points[x,0] = (x + 1) * midlen
	for x in range(19999):
		s1, i1, r1, p1, se1 = stats.linregress(x_data[:init_points[x,0]], y_data[:init_points[x,0]])
		s2, i2, r2, p2, se2 = stats.linregress(x_data[init_points[x,0]:], y_data[init_points[x,0]:])
		init_points[x,1] = r1**2
		init_points[x,2] = r2**2
		init_points[x,3] = (r1**2 + r2**2)
		print x, "completed"
	maxPos = np.argmax(init_points[900:19000,3])
	if init_points[maxPos-1,3] > init_points[maxPos+1,3]:
		bound = maxPos-1
	else:
		bound = maxPos+1
	plt.figure()
	plt.plot(init_points[:,3])
	plt.show()

#filePath = "C:\Users\ekram\Desktop\Spec4-0006-output.txt"
#filePath = "C:\Users\ekram\Desktop\Spec4-0110-output.txt"
filePath = R"C:\Users\ekram\Desktop\testFolder\Spec4-0474-output.txt"

## generate and break-up data

dataFile = np.genfromtxt(filePath, skip_header=1)

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
        if abs(int_points[x]) == 0:
            bound_pts[bndCnt,0] = x_points[x]
            bound_pts[bndCnt,1] = int_points[x]
            slopeSwitch = False
            bndCnt = bndCnt + 1
    else:
        if abs(int_points[x]) > 0:
            bound_pts[bndCnt,0] = x_points[x-1]
            bound_pts[bndCnt,1] = int_points[x-1]
            slopeSwitch = True
            bndCnt = bndCnt + 1
    if bndCnt > 4:
        break
print "bound_pts = "
print bound_pts

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
        bndCnt=bndCnt+1
    if bndCnt > 4:
        break

## Rescaled vectors to simplify following functions
approachT = timeCh1[Vbounds[1,2]:Vbounds[2,2]]
approachZ = 4.77*13*distance[Vbounds[1,2]:Vbounds[2,2]]
approachD = deflection[Vbounds[1,2]:Vbounds[2,2]]
approachF = k_L * approachD
retractT = timeCh1[Vbounds[3,2]:Vbounds[4,2]]
retractZ = 4.77*13*distance[Vbounds[3,2]:Vbounds[4,2]]
retractD = deflection[Vbounds[3,2]:Vbounds[4,2]]
retractF = k_L * retractD

multiLinReg(retractZ,retractF)

## Derivative of Approach curve, with rescaled Z-position
numPoints1 = len(approachZ)
points1 = 200
x1_points = np.zeros(points1)
y1_points = np.zeros(points1)
int_points1 = np.zeros(points1)
Abound_pts = np.zeros((4,3))
for x1 in range(points1):
    low = x1*numPoints1/points1
    high = low + numPoints1/points1
    if high > numPoints1:
        high = numPoints1
    s1, i1, r1, p1, e1 = stats.linregress(
        approachZ[low:high], approachD[low:high])
    x1_points[x1] = np.mean(approachZ[low:high])
    y1_points[x1] = s1
    int_points1[x1] = round(s1)

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
for x2 in range(points2):
    low = x2*numPoints2/points2
    high = low + numPoints2/points2
    if high > numPoints2:
        high = numPoints2
    s2, i2, r2, p2, e2 = stats.linregress(
        retractZ[low:high], retractD[low:high])
    x2_points[x2] = np.mean(retractZ[low:high])
    y2_points[x2] = s2
    int_points2[x2] = round(s2)

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
    if retractZ[x-1] > Rbound_pts[rBndCnt,1]:
        Rbounds[rBndCnt,0] = x
        Rbounds[rBndCnt,1] = retractZ[x]
        Rbounds[rBndCnt,2] = retractD[x]
        rBndCnt = rBndCnt + 1
    if rBndCnt > 3:
        break

# Find baseline & contact region linear regression
baselineS, baselineI, baselineR, baselineP, baselineE = stats.linregress(
        retractZ[Rbounds[1,0]:Rbounds[0,0]], retractD[Rbounds[1,0]:Rbounds[0,0]])
contactS, contactI, contactR, contactP, contactE = stats.linregress(
        retractZ[Rbounds[3,0]:Rbounds[2,0]], retractD[Rbounds[3,0]:Rbounds[2,0]])

x_shift = (baselineI - contactI)/(contactS - baselineS)
y_shift = contactS * x_shift + contactI

## Linear Regression on approach/retract regions
## __1 = slope ; __2 = intercept ; __3 = r_value ; __4 = p_value ; __5 = std_error
##
## Setup Speed
setup1,setup2,setup3,setup4,setup5 = stats.linregress(
    timeCh1[0:Vbounds[0,2]], distance[0:Vbounds[0,2]])
print "setup v =", abs(setup1), "nm/s"

## Approach Speed
appr1,appr2,appr3,appr4,appr5 = stats.linregress(
    approachT, approachZ)
print "approch v =", abs(appr1), "nm/s"

## Retract Speed
retr1,retr2,retr3,retr4,retr5 = stats.linregress(
    retractT, retractZ)
print "retract v =", abs(retr1), "nm/s"

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
plt.plot(retractZ,retractD)
plt.plot(retractZ - x_shift,retractD - y_shift)
plt.plot(0, 0, 'ro')
plt.plot(Rbounds[:,1], Rbounds[:,2], 'ro')
plt.ylabel("Deflection (nm)")
plt.xlabel("Z-position (V)") 
plt.axis([-100,100,-50,50])

plt.subplot(2,3,5)
plt.title("Approach")
plt.plot(approachZ,approachD)
plt.ylabel("Deflection (nm)")
plt.xlabel("Z-position (V)")

plt.subplot(2,3,4)
plt.title("Full Retract")
plt.plot(retractZ,retractD)
plt.plot(retractZ - x_shift,retractD - y_shift)
plt.plot(0, 0, 'ro')
plt.plot(Rbounds[:,1], Rbounds[:,2], 'ro')
plt.ylabel("Deflection (nm)")
plt.xlabel("Z-position (V)") 

plt.show()
plt.close()

