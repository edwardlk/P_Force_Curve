import time
from os import path, listdir, makedirs
from tkinter import Tk, filedialog
# import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp


def outputFiles(dataFiles, addon):
    L = []
    for x in range(len(dataFiles)):
        temp = dataFiles[x]
        L.append(temp[:-8] + addon)
    return L


def dataAlign(x, srcDir, dstDir, dataFiles, dataOutput):
    currentfile = dataFiles[x]
    outputfile = dataOutput[x]

    cFilePath = str(path.join(srcDir, currentfile))

    with open(str(cFilePath), 'r') as f:
        for i, l in enumerate(f):
            pass
        totalLines = i + 1

    foot1 = int(totalLines/2) - 3
    head2 = int(totalLines/2) + 11

    Ch1Table = np.genfromtxt(cFilePath, skip_header=11, skip_footer=foot1)
    DeflectTable = np.genfromtxt(cFilePath, skip_header=head2)

    combinedFile = np.hstack((Ch1Table, DeflectTable))

    np.savetxt(path.join(dstDir, outputfile), combinedFile,
               header="x(s) y0(V) x(s) y1(m)", comments="")
    print("Completed ", x+1, " of ", len(dataFiles), " files.")


def main():
    # Designate input and output directories.
    root = Tk()
    root.withdraw()

    info = ("Please select the folder that contains "
            "the data files you wish to reformat.")

    srcDir = filedialog.askdirectory(parent=root, initialdir="/", title=info)
    dstDir = path.join(srcDir, 'output')
    # csvDir = path.join(srcDir, 'csv')

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
    # if not path.exists(csvDir):
    #     makedirs(csvDir)
    # else:
    #     dataFiles.remove('csv')

    # Create output files' names
    dataOutput = outputFiles(dataFiles, '-output.txt')
    # csvOutput = outputFiles(dataFiles, '.csv')

    pool = mp.Pool(processes=5)
    for x in range(len(dataFiles)):
        pool.apply_async(dataAlign, args=(x, srcDir, dstDir, dataFiles,
                                          dataOutput,))
    pool.close()
    pool.join()
    print("Finished analyzing", path.split(srcDir)[1])
    print('It took {:.2f} seconds to analyze %d files.'.format(
          time.time()-start) % (len(dataFiles)))


if __name__ == '__main__':
    main()

# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
#
# x = DeflectTable[:, 0]
# y = DeflectTable[:, 1]*1000000000
#
# ax.scatter(x, y)
# plt.show()
# plt.close()
#
# # For grabing headers
# line_number = int(raw_input('Enter the line number: '))
# with open('packages.txt') as f:
#     i = 1  # might need i==1 for multiple lines
#     for line in f:
#         if i == line_number:
#             break
#         i += 1
#     # line now holds the line
#     # (or is empty if the file is smaller than that number)
#     print line
#
# for line in specFile:
#     print line
# specFile.close()
