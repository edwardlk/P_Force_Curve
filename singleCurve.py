from os import *
import math
import numpy as np
import matplotlib.pyplot as plt
# import Tkinter, tkFileDialog

# Import Data
# Col 1 - x-axis data in meters
# Col 2+ - alternating approach/retract curves, in nanometers
data = np.genfromtxt("C:\Users\Ed\Documents\GitHub\P_Force_Curve\Spec1.txt", \
skip_header = 20)

print "rows = ", data.shape[0]
print "columns = ", data.shape[1]
