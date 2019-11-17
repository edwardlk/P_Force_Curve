import time
from tkinter import Tk, filedialog
import multiprocessing as mp
from os import path, listdir, makedirs
from PFCfuncs import outputFiles, mainAnalysis


def main():
    # Designate input and output directories.

    testFile = False
    testMulti = True

    if testFile:
        srcDir = R"F:\TEST\fullAnalysisTest"
    else:
        root = Tk()
        root.withdraw()
        info = ("Please select the folder that contains the data files "
                "you wish to analyze.")
        srcDir = filedialog.askdirectory(parent=root,
                                         initialdir="/", title=info)
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
    # fillData = np.array([np.arange(len(dataFiles))]*7).T
    # df = pd.DataFrame(fillData, columns=['file', 'rupture force (pN)',
    #                                      'location', 'speed (nm/s)', 'WLC-P',
    #                                      'WLC-L0', 'x_off'])
    # df.to_pickle(path.join(csvDir, "dummy.pkl"))

    if testMulti:
        for x1 in range(len(dataFiles)):
            mainAnalysis(x1, srcDir, dstDir, csvDir, dataFiles, dataImg,
                         csvOutput, csvRupture)
    else:
        pool = mp.Pool(processes=5)
        for x1 in range(len(dataFiles)):
            pool.apply_async(mainAnalysis, args=(x1, srcDir, dstDir, csvDir,
                                                 dataFiles, dataImg, csvOutput,
                                                 csvRupture,))
        pool.close()
        pool.join()

    # df = pd.read_pickle(path.join(csvDir, "dummy.pkl"))
    # df.to_excel(path.join(csvDir, 'dataframe.xlsx'), sheet_name='Sheet1')
    print("Finished analyzing", path.split(srcDir)[1])
    print('It took {:.2f} seconds to analyze %d files.'.format(
          time.time()-start) % (len(dataFiles)))


if __name__ == '__main__':
    main()
