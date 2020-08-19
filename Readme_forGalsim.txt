Scripts to make GalSim chip files
(*** UNDER CONSTRUCTION ***)

First, run test_run.py.

---

Then make a FITS file with the noise information, e.g.:

python noise_run.py -i /fs/scratch/PCON0003/cond0007/SCA20829-noise/20191005_95k_1p1m0p1_noise_20829_001.fits -o temp2.fits -f 4 -n 100

The "-i" indicates the 1st input noise file
The "-o" indicates the 1st output file
The "-f" indicates the format number (4 for Roman data from the DCL)
The "-n" indicates the number of noise files to use (they are numbered sequentially)

---

Then run python assemble_galsimchipfile.py <script file>

The script file should be something like:
(comments can be added in lines starting with #)

# Label for the file
LABEL: Sample
# the SCA number
SCA: 20829
# input summary file from test_run.py
IN: /users/PCON0003/cond0088/Projects/detectors/sw_outputs/chris_20829st_summary.txt
# output file
OUT: out-chip/test.fits
# noise information
NOISE: temp2.fits

---

Right now, this makes the following HDUs:

primary
SOURCES     (information on the input files that went into this chip file
READ        (read noise HDU)
GAIN        (gain map HDU)

