import time
from tkinter import Tk, filedialog
import multiprocessing as mp
from os import path, listdir, makedirs
from PFCfuncs import outputFiles, mainAnalysis
import numpy as np
import pandas as pd


def main():
    # Designate input and output directories.

    testFile = False
    testMulti = True

    if testFile:
        # srcDir = R"F:\TEST\fullAnalysisTest"
        # srcDir = R'D:\TEST'
        srcDir = R'F:\TEST\errors'
        srcDir = R'D:\mlr_fix_test'
    else:
        root = Tk()
        root.withdraw()
        info = ("Please select the folder that contains the data files "
                "you wish to analyze.")
        srcDir = filedialog.askdirectory(parent=root,
                                         initialdir="/", title=info)
    try:
        k_L = float(input('Enter the cantilever stiffness (N/m): '))
        print('OK! Setting k_L = {} N/m.\nStarting analysis now...'.format(k_L))
    except Exception:
        k_L = 1.0
        print('Something went wrong, set k_L = 1.0 N/m.')
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
    outputPkl = path.join(csvDir, "dummy.pkl")
    outputCSV = path.join(csvDir, 'dataframe.xlsx')
    col_list = ['fnum', 'filename', 'min_location', 'fit_start', 'cStart', 'bStart',
                'Vb1', 'Vb2', 'Vb3', 'Vb4', 'Vb5']
    df = pd.DataFrame(columns=col_list)
    df.to_pickle(outputPkl)
    # df.to_excel(outputCSV, sheet_name='Sheet1')
    del df

    if testMulti:
        for x1 in range(len(dataFiles)):
            mainAnalysis(x1, k_L, srcDir, dstDir, csvDir, dataFiles, dataImg,
                         csvOutput, csvRupture, outputPkl)
    else:
        pool = mp.Pool(processes=5)
        for x1 in range(len(dataFiles)):
            pool.apply_async(mainAnalysis, args=(x1, k_L, srcDir, dstDir, csvDir,
                                                 dataFiles, dataImg, csvOutput,
                                                 csvRupture, outputPkl, ))
        pool.close()
        pool.join()

    outputDF = pd.read_pickle(outputPkl)
    outputDF.sort_values(by=['fnum'])
    outputDF.to_excel(outputCSV, sheet_name='Sheet1')

    print("Finished analyzing", path.split(srcDir)[1])
    print('It took {:.2f} seconds to analyze {} files.'.format(
          time.time()-start, len(dataFiles)))


if __name__ == '__main__':
    main()
