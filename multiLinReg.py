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

    # # Plots for testing
    # plt.figure()
    # plt.subplot(131)
    # plt.plot(init_points[:,1], 'r.')
    # plt.plot(init_points[:,2], 'b.')
    # plt.plot(init_points[:,3], 'y.')
    #
    # plt.subplot(132)
    # plt.plot(init_points[:,4], 'r.')
    # plt.plot(init_points[:,5], 'b.')
    # plt.plot(init_points[:,4] - init_points[:,5], 'b.')
    #
    # plt.subplot (133)
    # plt.plot(x_data,y_data)
    # plt.plot(x_data,s1*x_data + i1)
    # plt.plot(x_data,s2*x_data + i2)
    # plt.plot(x_data[init_points[maxPos, 0]],
    #          y_data[init_points[maxPos, 0]], 'ro')
    # plt.show()
    return s1, i1, s2, i2
