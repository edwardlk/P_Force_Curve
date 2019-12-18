import time
from os import path, listdir, makedirs
from tkinter import Tk, filedialog
import numpy as np
import multiprocessing as mp
import pandas as pd


def outputFiles(dataFiles, num, addon):
    ''' Adjust list of file names, removing 'num' characters and adding
        'addon'.
        Returns list.
    '''
    L = []
    for x in range(len(dataFiles)):
        temp = dataFiles[x]
        L.append(temp[:-num] + addon)
    return L


def dataAlign(x, srcDir, dstDir, dataFiles, dataOutput):
    ''' Data file realignment procedure for use in multiprocessing
        module.
    '''
    currentfile = dataFiles[x]
    outputfile = dataOutput[x]

    cFilePath = str(path.join(srcDir, currentfile))

    with open(str(cFilePath), 'r') as f:
        for i, l in enumerate(f):
            pass
        totalLines = i + 1

    foot1 = int(totalLines/2) - 4
    head2 = int(totalLines/2) + 11

    Ch1Table = np.genfromtxt(cFilePath, skip_header=11, skip_footer=foot1)
    DeflectTable = np.genfromtxt(cFilePath, skip_header=head2)

    combinedFile = np.hstack((Ch1Table, DeflectTable))

    ch1ys = ''
    deflTys = ''
    for x1 in range(Ch1Table.shape[1] - 1):
        ch1ys = ch1ys + ' y{}(V)'.format(x1)
        deflTys = deflTys + ' y{}(nm)'.format(x1)

    headtxt ='x(m)' + deflTys

    np.savetxt(path.join(dstDir, outputfile), DeflectTable,
               header=headtxt, comments="")
    print("Completed ", x+1, " of ", len(dataFiles), " files.")


def main():
    # Designate input and output directories.
    root = Tk()
    root.withdraw()

    info = ("Please select the folder that contains "
            "the data files you wish to reformat.")

    srcDir = filedialog.askdirectory(parent=root, initialdir="/", title=info)
    dstDir = path.join(srcDir, 'forcecurves')

    start = time.time()

    # Get file list from source directory
    dataFiles = listdir(srcDir)
    dataFiles.sort()

    # Make output directory if it does not exist
    # if the directory does exist, deletes 'output' from dataFiles
    if not path.exists(dstDir):
        makedirs(dstDir)
    else:
        dataFiles.remove('output')

    # Create output files' names
    dataOutput = outputFiles(dataFiles, 8, '-output.txt')

    # Convert files - multiprocessing
    pool = mp.Pool(processes=5)
    for x in range(len(dataFiles)):
        pool.apply_async(dataAlign, args=(x, srcDir, dstDir, dataFiles,
                                          dataOutput,))
    pool.close()
    pool.join()

    print("Finished analyzing", path.split(srcDir)[1])
    print('It took {:.1f} seconds to analyze {} files.'.format(
          time.time()-start, len(dataFiles)))


if __name__ == '__main__':
    main()
