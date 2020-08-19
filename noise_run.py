import sys
import re
import numpy
import time
from scipy import linalg
from astropy.io import fits
from os import path

import pyirc

# Get command-line arguments
fileformat = 1
nfile = 2
infile = outfile = ''
t_init = 2
for i in range(1, len(sys.argv)-1):
  if sys.argv[i] == '-f': fileformat = int(sys.argv[i+1])
  if sys.argv[i] == '-i': infile = sys.argv[i+1]
  if sys.argv[i] == '-o': outfile = sys.argv[i+1]
  if sys.argv[i] == '-n': nfile = int(sys.argv[i+1])
  if sys.argv[i] == '-t': t_init = int(sys.argv[i+1])

if infile == '':
  print('Error: no infile. Use -i option.')
  exit()

if outfile == '':
  print('Error: no outfile. Use -o option.')
  exit()

print('Input:', infile, ' +', nfile)
print('  ... Format =', fileformat)
print('Output:', outfile)

# Get array of new files:
m = re.search(r'^(.+)/(\d\d\d\d\d\d\d\d)([^/]+)\_(\d+).fits', infile)
if m:
  prefix = m.group(1)
  stamp = int(m.group(2))
  body = m.group(3)
else:
  print('Match failed.')
  exit()
filelist = []
for k in range(nfile):
  thisname = prefix + '/{:d}'.format(stamp) + body + '_{:03d}.fits'.format(k+1)
  count=0
  while not path.exists(thisname):
    count+=1
    thisname = prefix + '/{:d}'.format(stamp+count) + body + '_{:03d}.fits'.format(k+1)
    if count==10000:
      print('Failed to find file', k)
      exit()
  filelist.append(thisname)

# Now get the median of the first stable frames
nside = pyirc.get_nside(fileformat)
firstframe = numpy.zeros((nfile, nside, nside), dtype=numpy.uint16)
print(numpy.shape(firstframe))
for k in range(nfile):
  firstframe[k,:,:] = pyirc.load_segment(filelist[k], fileformat, [0,nside,0,nside], [t_init], False)

# Median and IQR
medarray = numpy.percentile(firstframe, 50, axis=0, interpolation='linear')
iqrarray = numpy.percentile(firstframe, 75, axis=0, interpolation='linear') - numpy.percentile(firstframe, 25, axis=0, interpolation='linear')

del firstframe

cdsnoise = numpy.zeros((nside,nside))
nband = 32
ntslice = pyirc.get_num_slices(fileformat, filelist[k])
nstep = ntslice//2 - 2
print('nstep =', nstep)
for j in range(nband):
  bandData = numpy.zeros((nfile, nstep, nside//nband, nside), dtype=numpy.int32)
  ymin = j*nside//nband
  ymax = ymin + nside//nband
  for k in range(nfile):
    print('band', j, 'file', k); sys.stdout.flush()
    for kb in range(nstep):
      bandData[k,kb,:,:] = pyirc.load_segment(filelist[k], fileformat, [0,nside,ymin,ymax], [t_init+2*kb  ], False).astype(numpy.int32)\
                             - pyirc.load_segment(filelist[k], fileformat, [0,nside,ymin,ymax], [t_init+2*kb+1], False).astype(numpy.int32)
  cdsnoise[ymin:ymax,:] = numpy.percentile(bandData, 75, axis=(0,1), interpolation='linear')\
                            - numpy.percentile(bandData, 25, axis=(0,1), interpolation='linear')
  del bandData
cdsnoise /= 1.34896

# Get PCAs
# split into this many channels
nch = 32
#
cdsFL = numpy.zeros((nfile, nside, nside), dtype=numpy.float32)
print(numpy.shape(cdsFL))
for k in range(nfile):
  print('build CDS ', k); sys.stdout.flush()
  cdsFL[k,:,:] = pyirc.load_segment(filelist[k], fileformat, [0,nside,0,nside], [t_init], False).astype(numpy.float32)\
                 - pyirc.load_segment(filelist[k], fileformat, [0,nside,0,nside], [ntslice-2], False).astype(numpy.float32)
  # reference pixel removal -- first sides, then top
  for j in range(nside):
    cdsFL[k,j,:] -= numpy.median(numpy.concatenate((cdsFL[k,j,:4],cdsFL[k,j,-4:])))
  for i in range(nch):
    xmin = i*nside//nch
    xmax = xmin + nside//nch
    cdsFL[k,:,xmin:xmax] -= numpy.median(cdsFL[k,-4:,xmin:xmax])
med_cdsFL = numpy.percentile(cdsFL, 50, axis=0, interpolation='linear').astype(numpy.float32)
iqr_cdsFL = (numpy.percentile(cdsFL, 75, axis=0, interpolation='linear') - numpy.percentile(cdsFL, 25, axis=0, interpolation='linear')).astype(numpy.float32)
iqr_cdsFL_med = numpy.median(iqr_cdsFL)
niqr_clip = 3.0
for k in range(nfile):
  cdsFL[k,:,:] -= med_cdsFL
  cdsFL[k,:,:] = numpy.where(cdsFL[k,:,:]>niqr_clip*iqr_cdsFL_med, niqr_clip*iqr_cdsFL_med, cdsFL[k,:,:])
  cdsFL[k,:,:] = numpy.where(cdsFL[k,:,:]<-niqr_clip*iqr_cdsFL_med, -niqr_clip*iqr_cdsFL_med, cdsFL[k,:,:])
#
# get covariance
C = numpy.zeros((nside, nside))
for i in range(nfile):
  print('cdsFL[i,:,:] =', cdsFL[i,:,:])
  for j in range(nfile):
    C[j,i] = numpy.sum(cdsFL[j,:,:].astype(numpy.float64)*cdsFL[i,:,:].astype(numpy.float64))
C /= float(nfile-1)
w,v = linalg.eigh(C)
max_eigenval = w[-1]
pca0_eigenvec = v[:,-1]
print('max eigenvalue =', max_eigenval)
print('eigenvec =', pca0_eigenvec)
print('normalization of eigenvec =', numpy.sum(pca0_eigenvec**2))
PCA0 = numpy.zeros((nside,nside))
for k in range(nfile):
  PCA0 += pca0_eigenvec[k] * cdsFL[k,:,:].astype(numpy.float64)
PCA0 /= numpy.sqrt(numpy.mean(PCA0**2))
#
pca0_stdev = numpy.sqrt(max_eigenval)/float(nside)
print('stdev of pca0 =', pca0_stdev)
# clean up memory
del med_cdsFL; del iqr_cdsFL; del cdsFL

# Generate output cube
outcube = numpy.zeros((4,nside,nside))
outcube[0,:,:] = 65535 - medarray
outcube[1,:,:] = iqrarray/1.34896
outcube[2,:,:] = cdsnoise/1.34896
outcube[3,:,:] = PCA0

primary_hdu = fits.PrimaryHDU()
out_hdu = fits.ImageHDU(outcube)
out_hdu.header['EXTNAME'] = 'NOISE'
out_hdu.header['PCA0_AMP'] = pca0_stdev
for k in range(nfile): out_hdu.header['NR_MF{:03d}'.format(k+1)] = filelist[k].strip()
out_hdu.header['NR_DATE'] = time.asctime(time.localtime(time.time()))

hdul = fits.HDUList([primary_hdu, out_hdu])
hdul.writeto(outfile, overwrite=True)

