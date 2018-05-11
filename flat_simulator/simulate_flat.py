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
"""
import sys
import numpy as np
from numpy.random import randn,poisson
import astropy.io.fits as fits
import fitsio
from fitsio import FITS,FITSHDR

# data cube attributes
nx, ny = 4096, 4096
#nx, ny = 128, 128
tsamp = 66
I = 2 # arbitrary scalar for now (e/s/pixel)
delta_t = 3 # arbitrary for now (s)
data_cube_Q = np.zeros((tsamp, nx, ny))
data_cube_S = np.zeros_like(data_cube_Q)
gain = 1.5 # arbitrary scalar e-/DN

# Start with 0 charge in the first frame (t=0)
for tdx in range(1, tsamp):
    mean = I*delta_t
    data_cube_Q[tdx,:,:] = data_cube_Q[tdx-1,:,:]+np.random.poisson(
        mean, data_cube_Q[tdx,:,:].shape)
    data_cube_S[tdx,:,:] =  data_cube_Q[tdx,:,:]/gain

# Open up an example DCL flat file and save the data cube
dclfile = 'Set_001_Test_0002.fits'
fitsio.write(dclfile, data_cube_S, clobber=True)

# Mean of a given slice checks out
# data_cube[1,:,:].mean()
# Try compression of data cube into file
# DCL file saved in 16-bit unsigned integers (look at header)
# End of script

"""
Things planned:
 * unsigned 16-bit integers
 * offset & clipping
 * I, dt -> floating point
 * reference pixels (4 around edge for all but WFC3 which has 5)
 * generate flats and darks
 * use real dark cube as read noise is reasonable, won't do hot pixels correctly, but ok for now

 * multiply I by quantum efficiency 

 * time step for read from t_a to t_a+1, need some configurable number of substeps so we can check for convergence
