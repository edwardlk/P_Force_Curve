import time
import datetime
from tkinter import Tk, filedialog
import multiprocessing as mp
from os import path, listdir, makedirs
from PFCfuncs import outputFiles, mainAnalysis, mainAnalysis2
import pandas as pd


def main():
    # Designate input and output directories.

    testFile = False
    testMulti = True

    if testFile:
        srcDir = R'F:\_data\Avidin-Biotin\2018-08-02_av-bioData\separated_files'
    else:
        root = Tk()
        root.withdraw()
        info = ("Please select the folder that contains the data files "
                "you wish to analyze.")
        srcDir = filedialog.askdirectory(parent=root,
                                         initialdir="/", title=info)
    try:
        k_L = float(input('Enter the cantilever stiffness (N/m): '))
        print('OK! Setting k_L = {} N/m.'
              '\nStarting analysis now...'.format(k_L))
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
    col_list = ['fnum', 'filename', 'min_location', 'fit_start', 'cStart',
                'bStart', 'v_r', 'Vb1', 'Vb2', 'Vb3', 'Vb4', 'Vb5']
    df = pd.DataFrame(columns=col_list)
    df.to_pickle(outputPkl)
    # df.to_excel(outputCSV, sheet_name='Sheet1')
    del df

    if testMulti:
        for x1 in range(len(dataFiles)):
            # mainAnalysis(x1, k_L, srcDir, dstDir, csvDir, dataFiles, dataImg,
            #              csvOutput, csvRupture, outputPkl)
            mainAnalysis2(x1, k_L, srcDir, dstDir, csvDir, dataFiles, dataImg,
                          csvOutput, csvRupture, outputPkl)
    else:
        pool = mp.Pool(processes=5)
        for x1 in range(len(dataFiles)):
            pool.apply_async(mainAnalysis, args=(x1, k_L, srcDir, dstDir,
                                                 csvDir, dataFiles, dataImg,
                                                 csvOutput, csvRupture,
                                                 outputPkl, ))
        pool.close()
        pool.join()

    outputDF = pd.read_pickle(outputPkl)
    outputDF.sort_values(by=['fnum'])
    outputDF.to_excel(outputCSV, sheet_name='Sheet1')

    totalTime = str(datetime.timedelta(seconds=int(time.time()-start)))

    print("Finished analyzing", path.split(srcDir)[1])
    print('It took {} to analyze {} files.'.format(totalTime, len(dataFiles)))


if __name__ == '__main__':
    main()
