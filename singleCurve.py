from os import *
import math
import numpy as np
import matplotlib.pyplot as plt
# import Tkinter, tkFileDialog

# Import Data
# Row 0 - x-axis data in meters, then convert to nm
# Row 1+ - alternating approach/retract curves, in nanometers
data = np.genfromtxt("C:\Users\Ed\Documents\GitHub\P_Force_Curve\Spec1.txt", \
                     skip_header = 20)
data_t = np.transpose(data)
data_t[0] = 10**9 * data_t[0]

Fcurves = (data_t.shape[0]-1)/2

print data_t[2].argmax()

plt.figure()
plt.plot(data_t[0], data_t[2], 'r')

plt.show()
plt.close()
