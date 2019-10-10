import time
from tkinter import Tk, filedialog
import numpy as np
import pandas as pd
import multiprocessing as mp
from os import path, listdir, makedirs
from PFCfuncs import outputFiles, fitAnalysis


def main():

    testing = True

    if testing:
        srcDir = R"D:/Ed/Desktop/New folder/fitTest"
        rupFileLoc = R"D:/Ed/Desktop/New folder/fitTest/proteinFit.csv"
    else:
        # Designate input and output directories.
        root = Tk()
        root.withdraw()

        info = ("Please select the folder that contains the data files "
            "you wish to analyze.")
        srcDir = filedialog.askdirectory(parent=root, initialdir="/", title=info)

        info2 = ("Please select the file that contains the locations of the "
                 "force rupture events you wish to analyze.")
        rupFileLoc = filedialog.askopenfilename(parent=root, initialdir=srcDir,
                                                title=info2)

    dstDir = path.join(srcDir, 'fitImages')
    csvDir = path.join(srcDir, 'fitCSVs')
    rupFile = path.split(rupFileLoc)[1]

    start = time.time()

    # Get file list from source directory
    dataFiles = listdir(srcDir)
    dataFiles.sort()

    # Remove non-datafiles from the data file list
    if rupFile in dataFiles:
        dataFiles.remove(rupFile)

    if not path.exists(dstDir):
        makedirs(dstDir)
    else:
        dataFiles.remove('fitImages')
    if not path.exists(csvDir):
        makedirs(csvDir)
    else:
        dataFiles.remove('fitCSVs')

    # Create output files' names
    rupImg = outputFiles(dataFiles, '-fit.png')
    rupOutput = outputFiles(dataFiles, '-fit.csv')

    rupGuess = pd.read_csv(rupFileLoc)

    print(rupGuess)

    # Pandas DataFrame to store measurements
    fillData = np.array([np.arange(len(rupGuess))]*6).T
    df = pd.DataFrame(fillData, columns=['file', 'rupture force (pN)',
                                         'location', 'WLC-P', 'WLC-L0', 'x_off'])
    df.to_pickle(path.join(csvDir, "dummy.pkl"))

    for x in range(len(rupGuess)):
        fitAnalysis(x, srcDir, dstDir, csvDir, rupGuess, dataFiles, rupImg, rupOutput)

    # pool = mp.Pool(processes=5)
    # for x in range(len(dataFiles)):
    #     pool.apply_async(fitAnalysis, args=(x, dataFiles,))
    # pool.close()
    # pool.join()

    df = pd.read_pickle(path.join(csvDir, "dummy.pkl"))
    df.to_excel(path.join(csvDir, 'dataframe.xlsx'), sheet_name='Sheet1')
    print("Finished analyzing", path.split(srcDir)[1])
    print('It took {:.2f} seconds to analyze %d files.'.format(
          time.time()-start) % (len(dataFiles)))


if __name__ == '__main__':
    main()