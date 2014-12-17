from os import *
import math
import numpy as np
import matplotlib.pyplot as plt
# import Tkinter, tkFileDialog

# Import Data
# Row 0 - x-axis data in meters, then convert to nm
# Row 1+ - alternating approach/retract curves, in nm
data = np.genfromtxt("C:\Users\Ed\Documents\GitHub\P_Force_Curve\data0417.txt", \
                     skip_header = 20)
data_t = np.transpose(data)
data_t[0] = 10**9 * data_t[0]

numCurves = (data_t.shape[0]-1)/2

maxPos = data_t[2].argmax()
minPos = data_t[2].argmin()

data_t[4] = data_t[4]/0.059
data_t[4] = data_t[4] - data_t[4][0]

Force = data_t[4] * 0.035

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel('Distance (nm)')
ax1.set_ylabel('Deflection (nm)', color='b')
ax1.plot(data_t[0], data_t[4], 'b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')

ax2 = ax1.twinx()
ax2.set_ylabel('Force (nN)', color='r')
ax2.plot(data_t[0], Force, 'r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')

plt.show()
plt.close()
