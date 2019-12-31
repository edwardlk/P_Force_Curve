import time
import datetime
from tkinter import Tk, filedialog
import pandas as pd
import multiprocessing as mp
from os import path, listdir, makedirs
from PFCfuncs import fitOutputFiles, fitAnalysis
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def distPlotting(csvDir, speedsList, df, modelName):
    numSpeeds = len(speedsList)

    if type(modelName) == str:
        mName = modelName.replace(' ', '')
    else:
        mName = modelName

    for x in range(numSpeeds):
        df1 = df[df['model#'] == modelName]
        df1 = df1[df1['v(nm/s)'] == speedsList[x]]
        sns.distplot(df1['rupture force (pN)'], kde=False)
        plt.xlabel('rupture force (pN)', fontsize=15)
        plt.ylabel("Frequency", fontsize=15)
        plt.title('Rupture Forces @ {} nm/s'.format(speedsList[x]))
        plt.savefig(path.join(csvDir, 'F_r-{}-{:.0f}nmps'.format(mName, speedsList[x])))
        plt.close()

    for x in range(numSpeeds):
        df1 = df[df['model#'] == modelName]
        df1 = df1[df1['v(nm/s)'] == speedsList[x]]
        sns.distplot(df1['L_C'], kde=False)
        plt.xlabel('L_C', fontsize=15)
        plt.ylabel("Frequency", fontsize=15)
        plt.title('Contour Lengths @ {} nm/s'.format(speedsList[x]))
        plt.savefig(path.join(csvDir, 'L_C-{}-{:.0f}nmps'.format(mName, speedsList[x])))
        plt.close()

    for x in range(numSpeeds):
        df1 = df[df['model#'] == modelName]
        df1 = df1[df1['v(nm/s)'] == speedsList[x]]
        sns.distplot(df1['L_P'], kde=False)
        plt.xlabel('L_P', fontsize=15)
        plt.ylabel("Frequency", fontsize=15)
        plt.title('Persistence Lengths @ {} nm/s'.format(speedsList[x]))
        plt.savefig(path.join(csvDir, 'L_P-{}-{:.0f}nmps'.format(mName, speedsList[x])))
        plt.close()


def main():

    testing = True
    testingMulti = True

    if testing:
        srcDir = R"F:\_data\Avidin-Biotin\2018-08-02_av-bioData\separated_files\RetractCSVs"
        rupFileLoc = (
            R"F:\_data\Avidin-Biotin\2018-08-02_av-bioData\separated_files\RetractCSVs\dataframe.xlsx")
    else:
        # Designate input and output directories.
        root = Tk()
        root.withdraw()

        info = ("Please select the folder that contains the data files "
                "you wish to analyze.")
        srcDir = filedialog.askdirectory(parent=root, initialdir="/",
                                         title=info)
        info2 = ("Please select the file that contains the locations of the "
                 "force rupture events you wish to analyze.")
        rupFileLoc = filedialog.askopenfilename(parent=root, initialdir=srcDir,
                                                title=info2)

    imgDir = path.join(srcDir, 'fitImages')
    csvDir = path.join(srcDir, 'fitCSVs')
    rupFile = path.split(rupFileLoc)[1]

    start = time.time()

    # Get file list from source directory
    dataFiles = listdir(srcDir)
    dataFiles.sort()

    # Remove non-datafiles from the data file list
    if rupFile in dataFiles:
        dataFiles.remove(rupFile)
    if 'dummy.pkl' in dataFiles:
        dataFiles.remove('dummy.pkl')

    if not path.exists(imgDir):
        makedirs(imgDir)
    else:
        dataFiles.remove('fitImages')
    if not path.exists(csvDir):
        makedirs(csvDir)
    else:
        dataFiles.remove('fitCSVs')

    rupGuess = pd.read_excel(rupFileLoc, index_col=0)

    # Create output files' names
    rupImg = fitOutputFiles(rupGuess, '-fit.png')
    rupOutput = fitOutputFiles(rupGuess, '-fit.csv')

    # Pandas DataFrame to store measurements
    col_list = ['file', 'model#', 'v(nm/s)', 'rupture force (pN)', 'location',
                'fit_method', 'function_evals', 'data_pts', 'variables',
                'chi-squared', 'r_chi-squared', 'Akaike_ic', 'Bayesian_ic',
                'L_P', 'L_P-err', 'L_C', 'L_C-err']
    # fillData = np.array([np.arange(len(rupGuess)*3)]*16).T
    fit_df = pd.DataFrame(columns=col_list)
    fit_df.to_pickle(path.join(csvDir, "dummy.pkl"))

    if testingMulti:
        for x in range(len(rupGuess)):
            fitAnalysis(x, srcDir, imgDir, csvDir, rupGuess, dataFiles, rupImg,
                        rupOutput)
    else:
        pool = mp.Pool(processes=5)
        for x in range(len(dataFiles)):
            pool.apply_async(fitAnalysis, args=(x, dataFiles,))
        pool.close()
        pool.join()

    df = pd.read_pickle(path.join(csvDir, "dummy.pkl"))
    df.to_excel(path.join(csvDir, 'dataframe.xlsx'), sheet_name='Sheet1')

    # plot histograms
    speedsList = df['v(nm/s)'].unique()
    modelList = df['model#'].unique()
    for x in range(len(modelList)):
        distPlotting(csvDir, speedsList, df, modelList[x])

    totalTime = str(datetime.timedelta(seconds=int(time.time()-start)))

    print("Finished analyzing", path.split(srcDir)[1])
    print('It took {} (H:M:S) seconds to analyze {} files.'.format(
        totalTime, len(dataFiles)))


if __name__ == '__main__':
    main()
