""" Start with a realization of a perfect detector
1.  The mean charge <Q_a(i,j)> = It_a
2.  Realization of a 4096 x 4096 pixel^2 grid with 66 time samples
This final result data cube will then be [4k, 4k, 66] in dimensions
and will require some kind of identifying header, currently taken
from one of the DCL flats.  This is charge.
3.  Convert to ADU/DN via gain ~ number(e-)/counts;  set ~ FWD/2^16

Run tests with 128x128 or something small...
current I is per pixel (units: e/s)
time t (units: s)
time step for read from t_a to t_a+1 (will need to check # for convergence)

NOTE: to run this, one needs: 
 ** a copy of a DCL flat file with the name set below
"""
import sys
import numpy as np
from numpy.random import randn,poisson
import astropy.io.fits as fits
import fitsio
from fitsio import FITS,FITSHDR
import sys
import re
sys.path.insert(0, '../')
#sys.path.insert(0, '/users/PCON0003/cond0080/src/solid-waffle/')
from pyirc import *

# Defaults
formatpars = 1
tsamp = 66
substep = 2
I = 10.0
QE = 0.8
delta_tsamp = 3.0 # arbitrary for now (s)
gain = 1.5 # arbitrary scalar e-/DN
outfile = 'DefaultOutput.fits'

# Read in information
config_file = sys.argv[1]
with open(config_file) as myf: content = myf.read().splitlines()
for line in content:
  # Format
  m = re.search(r'^FORMAT:\s*(\d+)', line)
  if m: formatpars = int(m.group(1))

  # Number of reads
  m = re.search(r'^NREADS:\s*(\d+)', line)
  if m: tsamp = int(m.group(1))
  # substeps
  m = re.search(r'^SUBSTEPS:\s*(\d+)', line)
  if m: substep = int(m.group(1))

  # Time step (s)
  m = re.search(r'^DT:\s*(\S+)', line)
  if m: delta_tsamp = float(m.group(1))

  # Gain (e/DN)
  m = re.search(r'^GAIN:\s*(\S+)', line)
  if m: gain = float(m.group(1))

  # Illumination (photons/s/pixel)
  m = re.search(r'^ILLUMINATION:\s*(\S+)', line)
  if m: I = float(m.group(1))
  # QE (Illumination * QE = current, e/s/pixel)
  m = re.search(r'^QE:\s*(\S+)', line)
  if m: QE = float(m.group(1))

  # Output file
  m = re.search(r'^OUTPUT:\s*(\S+)', line)
  if m: outfile = m.group(1)

# data cube attributes
N = nx = ny = get_nside(formatpars) # possibly redundant with nx,ny
# Reference pixels hard-coded to 4 rows/cols around border, true for
# all except WFC3 which has 5
xmin,xmax,ymin,ymax = 4,N-4,4,N-4 # Extent of non-reference pixels
nt_step = tsamp*substep # number of tot timesteps depending on convergence needs
delta_t = (delta_tsamp*tsamp)/nt_step # time between timesteps
allQ = np.zeros((substep, nx, ny))
data_cube_Q = np.zeros((tsamp, nx, ny))
data_cube_S = np.zeros_like(data_cube_Q)
count = 1

print 'side length =',N
print 'samples:', tsamp, 'x', delta_tsamp, 's; # substep =', substep
print 'Illumination:', I, 'ph/s/pix; QE =', QE

# Start with 0 charge in the first frame (t=0)
for tdx in range(1, nt_step):
  # Use either the flat current or a dark current (this setting will need
  # to be put into the config
  if light:
    mean = I*delta_t
  else:
    mean = I_dark*delta_t
    # Create realization of charge
    # This version uses less memory, but probably still sub-optimal
    idx = tdx%substep
    # Charge accumulates, dependent on quantum efficiency of the pixel
    allQ[idx,xmin:xmax,ymin:ymax] = allQ[idx-1,xmin:xmax,xmin:xmax] \
        +QE * np.random.poisson(mean, allQ[idx,xmin:xmax,xmin:xmax].shape)
    if (idx==0):
        data_cube_Q[count,:,:] = allQ[idx,:,:]
        allQ = np.zeros((substep, nx, ny))
        allQ[0,:,:] = data_cube_Q[count,:,:]
        count += 1

# Read in the read noise from a fits file generated with Bernie's ngxhrg
# Currently using one with one realization because the full one takes
# a long time to create
noisefile = 'ex_2.2.1.fits'
noise = fitsio.read(noisefile)
data_cube_Q[-1,:,:] += noise  # Adding to only the final time

# Convert charge to signal
data_cube_S = np.array(data_cube_Q/gain, dtype=np.uint16)

# Open up an example DCL flat file and save the data cube
#dclfile = 'Set_001_Test_0002.fits'
fitsio.write(outfile, data_cube_S, clobber=True)


# Mean of a given slice checks out
# data_cube[1,:,:].mean()
# Try compression of data cube into file
# DCL file saved in 16-bit unsigned integers (look at header)
# End of script

"""
Things planned:
 * offset & clipping
 * use real dark cube as read noise is reasonable, won't do hot pixels correctly, but ok for now
"""
