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

**TO DO:** Convert to python3
