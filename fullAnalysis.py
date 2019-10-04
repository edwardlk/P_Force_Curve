import sys
import time
from tkinter import Tk, filedialog
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import multiprocessing as mp
import multiLinReg as MLR
from os import path, listdir, makedirs
from scipy import stats
from scipy.optimize import curve_fit
from lmfit import Model
from polymerModels import WLCmodel, FJCmodel
from PFCfuncs import outputFiles, smooth, returnBoundaries, plotEverything, mainAnalysis


def main():
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
    df.to_excel(path.join(csvDir, 'dataframe.xlsx'), sheet_name='Sheet1')

    pool = mp.Pool(processes=5)
    for x1 in range(len(dataFiles)):
        pool.apply_async(mainAnalysis, args=(x1, srcDir, dstDir, csvDir,
                                             dataFiles, dataImg, csvOutput,
                                             csvRupture,))
    pool.close()
    pool.join()

    print("Finished analyzing", path.split(srcDir)[1])
    print('It took {:.2f} seconds to analyze %d files.'.format(time.time()-start)
          % (len(dataFiles)))


if __name__ == '__main__':
    main()
