import numpy as np
from scipy import stats


def multiLinReg(x_data, y_data):
    """ Inputs are ndarrays, x_data should be in nm's
    returns contact slope, contact intercept,
    baseline slope, baseline intercept"""

    # print "Running MLR"
    sepX = 1
    deltaX = 0.0
    while deltaX < 50:
        deltaX = abs(x_data[0] - x_data[sepX])
        sepX += 1
    midlen = len(x_data)/2000
    init_points = np.zeros((1999, 6))
    for x in range(1999):
        init_points[x, 0] = (x + 1) * midlen
    for x in range(1999):
        s1, i1, r1, p1, se1 = stats.linregress(x_data[:int(init_points[x, 0])],
                                               y_data[:int(init_points[x, 0])])
        s2, i2, r2, p2, se2 = stats.linregress(x_data[int(init_points[x, 0]):],
                                               y_data[int(init_points[x, 0]):])
        init_points[x, 1] = r1**2
        init_points[x, 2] = r2**2
        init_points[x, 3] = (r1**2 + r2**2)
        init_points[x, 4] = s1
        init_points[x, 5] = s2
        # if x % 100 == 0:
        #     print x, "completed"
    positons = init_points[:, 4] - init_points[:, 5]
    maxPos = int(np.argmax(positons[20:1980]) + 20)

    upperBound = int(init_points[maxPos, 0] - sepX/5)
    lowerBound = int(init_points[maxPos, 0] + sepX)

    if lowerBound > len(x_data):
        lowerBound = init_points[maxPos, 0]
    if upperBound > len(x_data):
        upperBound = init_points[maxPos, 0]

    # New fit curves, 1 for contact curve, 2 for baseline
    s1, i1, r1, p1, se1 = stats.linregress(x_data[:upperBound],
                                           y_data[:upperBound])
    s2, i2, r2, p2, se2 = stats.linregress(x_data[lowerBound:],
                                           y_data[lowerBound:])

    return s1, i1, s2, i2


def multiLinReg2(x_data, y_data):
    """ Attempts to find the contact point of the AFM curve by fitting 2 linear
        regressions, one for the baseline and one for the contact line.
        Returns contact slope, intercept, baseline slope, intercept.
    """
    iteration = 0

    left = 0
    right = len(x_data) - 1
    mid = (right + left) // 2
    midLast = 0

    while iteration < 100:
        if midLast == mid:
            break
        pt1 = (left + mid) // 2
        pt2 = (right + mid) // 2
        s1c, i1c, r1c, p1c, se1c = stats.linregress(x_data[:pt1],
                                                    y_data[:pt1])
        s1b, i1b, r1b, p1b, se1b = stats.linregress(x_data[pt1:],
                                                    y_data[pt1:])
        s2c, i2c, r2c, p2c, se2c = stats.linregress(x_data[:pt2],
                                                    y_data[:pt2])
        s2b, i2b, r2b, p2b, se2b = stats.linregress(x_data[pt2:],
                                                    y_data[pt2:])
        s1 = s1c - s1b
        s2 = s2c - s2b

        if s1 > s2:
            right = mid - 1
            sout = [s1c, i1c, s1b, i1b]
            midLast = mid
        else:
            left = mid + 1
            sout = [s2c, i2c, s2b, i2b]
            midLast = mid
        mid = (right + left) // 2
        iteration += 1
    return sout[0], sout[1], sout[2], sout[3]


def multiLinRegDirect(x_data, y_data, cStart, bStart):
    """ Attempts to find the contact point of the AFM curve by fitting 2 linear
        regressions, one for the baseline and one for the contact line.
        Returns contact slope, intercept, baseline slope, intercept.
    """
    # Create list of points for linear regressions
    regPts = [0, (len(x_data)-1)]

    # Add other linear regression points to list
    x_data1 = abs(x_data - cStart)
    cStartInt = x_data1.argmin()
    regPts.append(cStartInt)
    x_data2 = abs(x_data - bStart)
    bStartInt = x_data2.argmin()
    regPts.append(bStartInt)
    regPts.sort()

    cS, cI, r1c, p1c, se1c = stats.linregress(x_data[regPts[0]:regPts[1]],
                                              y_data[regPts[0]:regPts[1]])
    bS, bI, r1c, p1c, se1c = stats.linregress(x_data[regPts[2]:regPts[3]],
                                              y_data[regPts[2]:regPts[3]])

    return cS, cI, bS, bI
