## Protein Force Curve Analysis
The files in this document are a variety of different python (2.7) scripts that
I used over the years to format and analyze protein force curve data files
obtained via an Agilent AFM & RHK SPM1000 controller.

Below I'll list what each of the files was used for and any further work that
needs to be done to make the scripts more useful.

### fixfiles.py
Used to make data captured in 'Time Spec' easier to load by simplifying the
.SM3spec files. After using the RHK provided tools to convert the .sm4 files to
text files (.SM3spec), you'll end up with something like this:

    **Spec4-001.SM3spec**
    [header]
    time(s) Z-piezo(V)
    t0   V0
    t1   V1
    ...

    [header]
    time(s) deflection(nm)
    t0   d0
    t1   d1
    ...

This script will batch convert these to text files, removing the headers and
realigning the columns to look like this:

    **Spec4-001-output.txt**
    time(s) Z-piezo(V) time(s) deflection (nm)
    t0   V0    t0    d0
    t1   V1    t1    d1
    ...

To use, put only the text files that you want to convert into a folder, run the
script, when the dialog box opens select the folder that contains your data
files and then wait until it finishes. The converted files will be in a new
folder called 'output' that was created within the folder that you selected.

*TO DO:*
- Convert to python3

### fullAnalysis.py
Used to do a batch analysis of protein force curves captured in Time Spec mode.
Assumes data files have a format similar to the **Spec4-001-output.txt**
example given above, with one curve in each file.

A rough outline of the analysis it performs:
1. Uses Z-piezo position to discard liftoff & hold data, keeping only approach
& retract data.
2. Fits the contact & far-from-surface portions of the retract data to find the
surface location and remove any y-offset.
3. Calculates retract speed.
4. Smooth's the force curves with various wind sizes
5. Finds the rupture force & location, then tries to fit the data leading up to
the rupture using the worm-like chain (WLC) model.
6. Saves all the rupture data & WLC parameters to a separate dataframe
7. Creates 2 output files, one containing rescaled retract curve data, the other
with only the data that it trird to fit to the WLC model.
8. Prints various figure demonstrating how the data was transformed & fit, to
making it easier to find the files that the program failed to properly analyze.

*TO DO:*
- Convert to python3
- ask for cantilever stiffness
- allow switching b/w Time Spec & regular Spec files
- multiple Time Spec curves in a single file.
- make step 1 more robust for handling different Time Spec setups (e.g. w/ & w/o
  hold on surface)

### singleCurve.py
Testing script for what was eventually added to the fullAnalysis script. Kept
for reference.

### multiLinReg.py & polymerModels.py
Called in the analysis scripts. multiLinReg for finding the x- and y-offsets of
a retract curve. polymerModels has the WLC fit function.

*updated 7/20/18*
