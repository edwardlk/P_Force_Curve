import pandas as pd
from os import path, listdir, makedirs


def outputFiles(fileName, num, addon, numSpecs):
    ''' Adjust list of file names, removing 'num' characters and adding
        'addon'.
        Returns list.
    '''

    L = []
    for x in range(numSpecs):
        temp = fileName[:-num]
        L.append(temp + '-c{}'.format(str(x).zfill(3)) + addon)
    return L


srcDir = R'F:\_data\Avidin-Biotin\2018-08-09_av-bioData'

outDir = path.join(srcDir, 'separated_files')

# Get file list from source directory
dataFiles = listdir(srcDir)
dataFiles.sort()

if not path.exists(outDir):
    makedirs(outDir)
else:
    dataFiles.remove('separated_files')

rows2skip = [0, 1, 2, 3, 4, 5, 6, 7, 9, 10]

for x in range(len(dataFiles)):
    data_file = pd.read_table(path.join(srcDir, dataFiles[x]),
                              skipinitialspace=True, skiprows=rows2skip)

    numSpecs = data_file.shape[1] - 1
    out_files = outputFiles(dataFiles[x], 8, '.txt', numSpecs)

    col_names = data_file.columns

    for x1 in range(len(out_files)):
        col2write = [col_names[0], col_names[x1+1]]
        data_file.to_csv(path.join(outDir, out_files[x1]), columns=col2write,
                         index=False)

    print('finished {} of {}.'.format(x, len(dataFiles)))
