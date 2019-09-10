# Script to convert 4D to 3D images

import sys
import numpy
import time
import re
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits

# prototype: '../SCA20026-run274/may02_HG_smallLED_10P000_SN20026_'

stem = sys.argv[1]
istart = int(sys.argv[2])
n = int(sys.argv[3])

subfirst = True

# Gain settings:
nxg = 128 # channel width to use for gain determination
nchan = 4096//nxg # number of channels

nlDeg = 5 # degree of nonlinear fit, charge(DN)

dxchop = 4 # number of columns to remove on each end
dychop = 1 # number of rows to remove at bottom

graw = numpy.zeros((nchan))
icpt = numpy.zeros((nchan))
init_offset = numpy.zeros((nchan, n))
final_offset = numpy.zeros((nchan, n))
beta_gain = numpy.zeros((nchan, n))
poly_rev = numpy.zeros((nchan, n, nlDeg-1))

bdrop = 2 # number of frames to drop at beginning -- for settling
edrop = 4 # number of frames to drop at end -- sometimes there are data writing problems?

# random number generator for de-quantization of medians
numpy.random.seed(64) # this setting makes the result reproducible

for i in range(n):

  # Dark and flat files
  lightfile = stem + ('{:03d}'.format(istart+i)) + '.fits'
  print 'light [', i, ']: ', lightfile
  myCube = fits.open(lightfile)
  if edrop>0:
    data_3D = myCube[1].data[0,bdrop:-edrop,dychop:,:]
  else:
    data_3D = myCube[1].data[0,bdrop:,dychop:,:]
  data_3D = data_3D.astype(float); data_3D *= -1.
  data_3D += numpy.random.uniform(low=-0.5, high=0.5, size=numpy.shape(data_3D))
  if i==0:
    # setup -- number of time slices, photon transfer curve mean & var
    nt = numpy.shape(data_3D)[0]
    print 'number of time steps:', nt
    PTCmean = numpy.zeros((nt,nchan))
    PTCvar = numpy.zeros((nt,nchan))
    med_channel = numpy.zeros((nt,nchan,n))

  first_sub = numpy.copy(data_3D)
  for t in range(nt):
    first_sub[t,:,:] -= data_3D[0,:,:]
  darkfile = re.sub('_SN', '_dark_SN', lightfile)
  print 'dark:', darkfile
  myCube = fits.open(darkfile)
  if edrop>0:
    data_3D += myCube[1].data[0,bdrop:-edrop,dychop:,:].astype(float)
  else:
    data_3D += myCube[1].data[0,bdrop:,dychop:,:].astype(float)

  if (i==0):
    outname = re.sub(r'\.fits', r'_diff_nocds.fits', lightfile)
    print '>>', outname
    hdu = fits.PrimaryHDU(data_3D[:,:,:])
    hdul = fits.HDUList([hdu])
    #hdul.writeto(outname, overwrite=True)
    first_sub_i0 = first_sub

  # Get non-linearity curves
  for ichan in range(nchan):
    xmin = nxg*ichan+dxchop
    xmax = xmin+nxg-2*dxchop
    nl_y = numpy.zeros((nt))
    for t in range(nt):
      nl_y[t] = numpy.median(data_3D[t,:,xmin:xmax])
    p = numpy.poly1d(numpy.polyfit(numpy.array(range(nt)), nl_y, 2))
    beta_gain[ichan,i] = -p.c[-3]/(p.c[-2]**2-4*p.c[-1]*p.c[-3]) # nonlinearity parameter beta*g in DN^-1
    # note annoying fact that polynomial coefficients are reversed

    # alternative solution -- fit charge(signal) instead of signal(charge)
    pp = numpy.poly1d(numpy.polyfit(nl_y, numpy.array(range(nt)), nlDeg))
    for k in range(2,nlDeg+1):
      poly_rev[ichan,i,k-2] = pp.c[-k-1]/pp.c[-2]
    # polynomial coefficients for function(signal) that linearizes but has f(S) = 0 + S + [...]S^2 + ...

  # subtract first frame
  if subfirst:
    # save offsets for first frame used
    for ichan in range(nchan):
      xmin = nxg*ichan+dxchop
      xmax = xmin+nxg-2*dxchop
      init_offset[ichan,i] = numpy.median(data_3D[0,:,xmin:xmax])
      final_offset[ichan,i] = numpy.median(data_3D[-1,:,xmin:xmax])
    for t in range(1,nt): data_3D[t,:,:] -= data_3D[0,:,:]
    data_3D[0,:,:] = 0

  if (i==0):
    outname = re.sub(r'\.fits', r'_diff.fits', lightfile)
    print '>>', outname
    hdu = fits.PrimaryHDU(data_3D[:,:,:])
    hdul = fits.HDUList([hdu])
    #hdul.writeto(outname, overwrite=True)
    first_sub_i0 = first_sub

  if (i==1):
    # get gain
    diff_image = (first_sub-first_sub_i0)/numpy.sqrt(2.)
    mean_image = (first_sub+first_sub_i0)/2.
    for ichan in range(nchan):
      xmin = nxg*ichan+dxchop
      xmax = xmin+nxg-2*dxchop
      for t in range(nt):
        PTCmean[t,ichan] = numpy.median(mean_image[t,:,xmin:xmax])
        PTCvar[t,ichan] = ((numpy.percentile(diff_image[t,:,xmin:xmax],75.) - numpy.percentile(diff_image[t,:,xmin:xmax],25.))/1.349)**2
      A = numpy.vstack([PTCmean[1:,ichan], numpy.ones(nt-1)]).T
      graw[ichan], icpt[ichan] = numpy.linalg.lstsq(A, PTCvar[1:,ichan])[0]
      graw[ichan] = 1./graw[ichan]

  # compute medians, etc. on the difference image
  for ichan in range(nchan):
    xmin = nxg*ichan+dxchop
    xmax = xmin+nxg-2*dxchop
    for t in range(nt):
      med_channel[t,ichan,i] = numpy.median(data_3D[t,:,xmin:xmax])

# Print PTC table
thisOut = open(stem+'_ptc.txt', 'w')
for ichan in range(nchan):
  thisOut.write('{:12.5E} {:12.5E}   '.format(icpt[ichan], graw[ichan]))
  for t in range(nt):
    thisOut.write(' {:12.5E} {:11.5E}'.format(PTCmean[t,ichan], PTCvar[t,ichan]))
  thisOut.write('\n')
thisOut.close()

# Print median table
# each column is a channel; each row is a time step
thisOut = open(stem+'_med.txt', 'w')
for t in range(nt):
  for ichan in range(nchan):
    thisOut.write(' {:12.5E}'.format(numpy.median(med_channel[t,ichan,:])))
  thisOut.write('\n')
thisOut.close()

# Print mean & std dev on mean table
# each column is a channel; each row is a time step
thisOut = open(stem+'_mean.txt', 'w')
for t in range(nt):
  for ichan in range(nchan):
    thisOut.write(' {:12.5E} {:11.5E}'.format(numpy.mean(med_channel[t,ichan,:]), numpy.std(med_channel[t,ichan,:], ddof=1)/numpy.sqrt(n) ))
  thisOut.write('\n')
thisOut.close()

# diagnostics on ordering
thisOut = open(stem+'_argorder.txt', 'w')
for t in range(nt):
  for ichan in range(nchan):
    argvec = numpy.argsort(med_channel[t,ichan,:]) + 1
    for k in range(n):
      thisOut.write(' {:2d}'.format(argvec[k]))
  thisOut.write('\n')
thisOut.close()
# values in channel 1
thisOut = open(stem+'_vals_chan1.txt', 'w')
for t in range(nt):
  for i in range(n):
      thisOut.write(' {:12.5E}'.format(med_channel[t,0,i]))
  thisOut.write('\n')
thisOut.close()

# initial offsets
thisOut = open(stem+'_offsets.txt', 'w')
for ichan in range(nchan):
  thisOut.write('{:9.2f} {:9.2f} {:12.5E}  '.format(numpy.mean(init_offset[ichan,:]), numpy.mean(final_offset[ichan,:]), numpy.mean(beta_gain[ichan,:])))
  for k in range(2,nlDeg+1):
      thisOut.write(' {:19.12E}'.format(numpy.mean(poly_rev[ichan,:,k-2])))
  thisOut.write('\n');
thisOut.close()
print 'offsets: ', numpy.mean(init_offset[0,:]), numpy.mean(final_offset[0,:])

thisOut = open('currents.txt', 'a')
thisOut.write(stem + ' ')
thisOut.close
