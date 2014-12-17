from os import *
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# import Tkinter, tkFileDialog

def linearF(x, m, b):
    return m * x + b

# Import Data
# Row 0 - x-axis data in meters, then convert to nm
# Row 1+ - alternating approach/retract curves, in nm
data = np.genfromtxt("C:\Users\Ed\Documents\GitHub\P_Force_Curve\data0417.txt", \
                     skip_header = 20)
data_t = np.transpose(data)
data_t[0] = 10**9 * data_t[0]

numCurves = (data_t.shape[0]-1)/2

maxPos = data_t[4].argmax()
minPos = data_t[4].argmin()

# 0.059 is slope of contact curve, need to add code to find that first.
data_t[4] = data_t[4]/0.059

#Finding & fitting the baseline
popt, pcov = curve_fit(linearF, data_t[0][0:(minPos-1000)], data_t[4][0:(minPos-1000)])

# Plot curve, y_l showing deflection and y_r showing force
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel('Distance (nm)')
ax1.set_ylabel('Deflection (nm)')
ax1.plot(data_t[0], data_t[4]-popt[1], 'r')
y1, y2 = ax1.get_ylim()

ax2 = ax1.twinx()
ax2.set_ylabel('Force (nN)')
ax2.set_ylim(y1*0.035, y2*0.035)

plt.show()
plt.close()
