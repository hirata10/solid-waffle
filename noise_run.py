import sys
import re
import numpy
import time
from scipy import linalg
from scipy.signal import convolve
from astropy.io import fits
from os import path

import pyirc

# Get command-line arguments
fileformat = 1
nfile = 2
infile = outfile = ''
t_init = 2
rowh = 7
cds_cut = 10.
tnoiseframe = 55
refout = False
for i in range(1, len(sys.argv)):
  if sys.argv[i] == '-f': fileformat = int(sys.argv[i+1])
  if sys.argv[i] == '-i': infile = sys.argv[i+1]
  if sys.argv[i] == '-o': outfile = sys.argv[i+1]
  if sys.argv[i] == '-n': nfile = int(sys.argv[i+1])
  if sys.argv[i] == '-t': t_init = int(sys.argv[i+1])
  if sys.argv[i] == '-cd': cds_cut = float(sys.argv[i+1])
  if sys.argv[i] == '-rh': rowh = int(sys.argv[i+1])
  if sys.argv[i] == '-tn': tnoiseframe = int(sys.argv[i+1])
  if sys.argv[i] == '-ro': refout = True # has reference output

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

# get group time
with fits.open(filelist[0]) as G:
  tgroup = float(G[0].header['TGROUP'])
  if tgroup<.01: tgroup = float(G[0].header['TFRAME'])

# Now get the median of the first stable frames
nside = pyirc.get_nside(fileformat)
ntslice = pyirc.get_num_slices(fileformat, filelist[k])
firstframe = numpy.zeros((nfile, nside, nside), dtype=numpy.uint16)
print(numpy.shape(firstframe))
for k in range(nfile):
  firstframe[k,:,:] = pyirc.load_segment(filelist[k], fileformat, [0,nside,0,nside], [t_init], False)

# Median and IQR
medarray = numpy.percentile(firstframe, 50, axis=0, interpolation='linear')
iqrarray = numpy.percentile(firstframe, 75, axis=0, interpolation='linear') - numpy.percentile(firstframe, 25, axis=0, interpolation='linear')

del firstframe

# Dark maps
darkcds = numpy.zeros((nfile, nside, nside), dtype=numpy.float32)
for k in range(nfile):
  darkcds[k,:,:] = (pyirc.load_segment(filelist[k], fileformat, [0,nside,0,nside], [t_init], False).astype(numpy.float32)\
                      -pyirc.load_segment(filelist[k], fileformat, [0,nside,0,nside], [t_init+1], False).astype(numpy.float32))\
                      + ( (k+.5)/nfile - .5)
dark1 = numpy.median(darkcds, axis=0)/tgroup
dark1err = (numpy.percentile(darkcds, 75, axis=0) - numpy.percentile(darkcds, 25, axis=0))/tgroup * 0.9291/numpy.sqrt(nfile-1)
# comment: the "sigma on the median" of a Gaussian distribution is 0.9291 IQR
# (where IQR is the inter-quartile range)
for k in range(nfile):
  darkcds[k,:,:] = (pyirc.load_segment(filelist[k], fileformat, [0,nside,0,nside], [t_init], False).astype(numpy.float32)\
                      -pyirc.load_segment(filelist[k], fileformat, [0,nside,0,nside], [ntslice], False).astype(numpy.float32))\
       	       	      +	( (k+.5)/nfile - .5)
dark2 = numpy.median(darkcds, axis=0)/tgroup/(ntslice-t_init)
dark2err = (numpy.percentile(darkcds, 75, axis=0) - numpy.percentile(darkcds, 25, axis=0))/tgroup/(ntslice-t_init) * 0.9291/numpy.sqrt(nfile-1)
del darkcds

cdsnoise = numpy.zeros((nside,nside))
nband = 32
nstep = ntslice//2 - 2
print('nstep =', nstep)
for j in range(nband):
  bandData = numpy.zeros((nfile, nstep, nside//nband, nside), dtype=numpy.float32)
  ymin = j*nside//nband
  ymax = ymin + nside//nband
  for k in range(nfile):
    print('CDS band', j, 'file', k); sys.stdout.flush()
    for kb in range(nstep):
      bandData[k,kb,:,:] = pyirc.load_segment(filelist[k], fileformat, [0,nside,ymin,ymax], [t_init+2*kb  ], False).astype(numpy.float32)\
                             - pyirc.load_segment(filelist[k], fileformat, [0,nside,ymin,ymax], [t_init+2*kb+1], False).astype(numpy.float32)\
                           + ( (k+.5)/nfile + (kb+.5)/nstep -1.) # uniformly distribute to avoid quantization
  cdsnoise[ymin:ymax,:] = numpy.percentile(bandData, 75, axis=(0,1), interpolation='linear')\
                            - numpy.percentile(bandData, 25, axis=(0,1), interpolation='linear')
  del bandData
cdsnoise /= 1.34896

# total noise
tnoise = numpy.zeros((nside,nside))
nband = 32
for j in range(nband):
  bandData = numpy.zeros((nfile, nside//nband, nside))
  ymin = j*nside//nband
  ymax = ymin + nside//nband
  for k in range(nfile):
    print('TN band', j, 'file', k); sys.stdout.flush()
    for kb in range(tnoiseframe):
      bandData[k,:,:] = bandData[k,:,:]\
                         + (pyirc.load_segment(filelist[k], fileformat, [0,nside,ymin,ymax], [t_init+kb], False)\
                           + ((k+.5)/nfile - .5) )* (kb-(tnoiseframe-1)/2.) # weight by how far from center
  tnoise[ymin:ymax,:] = numpy.percentile(bandData, 75, axis=0, interpolation='linear')\
                            - numpy.percentile(bandData, 25, axis=0, interpolation='linear')
  del bandData
tnoise /= tnoiseframe*(tnoiseframe+1.)/12. # normalize so -1/2 --> +1/2 DN (change of 1) corresponds to one unit
  # this is sum (kb-(tnoiseframe-1)/2.) [ (kb-(tnoiseframe-1)/2.) / (tnoiseframe-1) ]
tnoise /= 1.34896

# ACN noise information
nband = 32 # number of channels
nfile2 = min(nfile,16)
ntslice = pyirc.get_num_slices(fileformat, filelist[k])
nstep = ntslice//2 - 2
print('nstep =', nstep)
acnband = numpy.zeros((nband))
aclip = 3.*1.34896*numpy.median(cdsnoise)
for j in range(nband):
  print('ACN, band', j)
  bandData = numpy.zeros((nfile2, nstep, nside, nside//nband))
  xmin = j*nside//nband
  xmax = xmin + nside//nband
  for k in range(nfile2):
    for kb in range(nstep):
      bandData[k,kb,:,:] = pyirc.load_segment(filelist[k], fileformat, [xmin,xmax,0,nside], [t_init+2*kb  ], False).astype(numpy.float64)\
                             - pyirc.load_segment(filelist[k], fileformat, [xmin,xmax,0,nside], [t_init+2*kb+1], False).astype(numpy.float64)
  for kb in range(nstep):
    bref = numpy.median(bandData[:,kb,:,:], axis=0)
    for k in range(nfile2): bandData[k,kb,:,:] -= bref
  # clip
  bandData = numpy.where(bandData>aclip, aclip, bandData)
  bandData = numpy.where(bandData<-aclip, -aclip, bandData)

  # get power spectrum
  bandPower = numpy.absolute(numpy.fft.fft(bandData, axis=3, norm='ortho'))**2
  print('shape of bandPower:', numpy.shape(bandPower))
  bandPowerDiff = bandPower[:,:,:,64] - numpy.mean(bandPower[:,:,:,60:64], axis=3)
  print('shape of bandPowerDiff:', numpy.shape(bandPowerDiff))
  #for i in range(128): print('{:3d} {:12.6f}'.format(i, numpy.mean(bandPower[:,:,:,i])))
  acnband[j] = numpy.mean(bandPowerDiff)/128./2. # divide by length; then by 2 to get per sample
  print('band', j, ': ACN =', acnband[j], 'DN**2'); sys.stdout.flush()
  del bandData
acnnoise = numpy.median(acnband) * nfile2/(nfile2-1.)
if acnnoise<0:
  acnnoise=0
else:
  acnnoise = numpy.sqrt(acnnoise)
print('overall ACN', acnnoise)

# 1/f noise information
nband = 32 # number of channels
extra_power_spectra = 1
if refout: extra_power_spectra = 3   # how many additional power spectra we compute

ntslice = pyirc.get_num_slices(fileformat, filelist[k])
nstep = ntslice//2 - 2
print('nstep =', nstep)
acnband = numpy.zeros((nband))
aclip = 3.*1.34896*numpy.median(cdsnoise)
avebandData = numpy.zeros((nfile2, nstep, nside, nside//nband))
slength = 4096*(4096//nband+rowh)
PS = numpy.zeros((nband+extra_power_spectra,slength))
for j in range(nband+extra_power_spectra):
  print('1/f, band', j); sys.stdout.flush()
  bandData = numpy.zeros((nfile2, nstep, nside, nside//nband))
  if j<nband:
    xmin = j*nside//nband
    xmax = xmin + nside//nband
    for k in range(nfile2):
      for kb in range(nstep):
        bandData[k,kb,:,:] = pyirc.load_segment(filelist[k], fileformat, [xmin,xmax,0,nside], [t_init+2*kb  ], False).astype(numpy.float64)\
                               - pyirc.load_segment(filelist[k], fileformat, [xmin,xmax,0,nside], [t_init+2*kb+1], False).astype(numpy.float64)
    for kb in range(nstep):
      bref = numpy.median(bandData[:,kb,:,:], axis=0)
      for k in range(nfile2): bandData[k,kb,:,:] -= bref
    # clip
    bandData = numpy.where(bandData>aclip, aclip, bandData)
    bandData = numpy.where(bandData<-aclip, -aclip, bandData)

    # flip even channels
    if j%2==1: bandData = numpy.flip(bandData, axis=3)
    
    avebandData += bandData/nband
  elif j==nband:
    bandData[:,:,:,:] = avebandData
  else:
    # now we are working with the reference output.
    # get reference +/- average of the other channels
    xmin = nside
    xmax = xmin + nside//nband
    for k in range(nfile2):
      for kb in range(nstep):
        bandData[k,kb,:,:] = pyirc.load_segment(filelist[k], fileformat, [xmin,xmax,0,nside], [t_init+2*kb  ], False).astype(numpy.float64)\
                               - pyirc.load_segment(filelist[k], fileformat, [xmin,xmax,0,nside], [t_init+2*kb+1], False).astype(numpy.float64)
    for kb in range(nstep):
      bref = numpy.median(bandData[:,kb,:,:], axis=0)
      for k in range(nfile2): bandData[k,kb,:,:] -= bref
    # clip
    bandData = numpy.clip(bandData,-aclip,aclip)

    if j==nband+1:
      bandData[:,:,:,:] += avebandData
    else:
      bandData[:,:,:,:] -= avebandData

  # get power spectrum
  mySeries = numpy.zeros((nfile2, nstep, slength))
  for r in range(4096): mySeries[:,:,r*(4096//nband+rowh):r*(4096//nband+rowh)+4096//nband] = bandData[:,:,r,:]

  bandPower = numpy.mean(numpy.absolute(numpy.fft.fft(mySeries, axis=2, norm='ortho'))**2, axis=(0,1))
  bandPower *= (1.+rowh*nband/4096.)/slength/2. # divide by length; then by 2 to get per sample instead of CDS
  print('shape of bandPower:', numpy.shape(bandPower))
  PS[j,:] = bandPower

  del bandData; del bandPower; del mySeries

for i in range(1024):
  print('{:6d} {:11.5E} {:12.5E}'.format(i, PS[nband,i], numpy.median(PS[:nband,i])-PS[nband,i]))

ncut = 4096
c_pink = numpy.sum(PS[nband,1:ncut+1]) - ncut*numpy.median(PS[nband,:])
u_pink = numpy.sum(PS[:nband,1:ncut+1])/nband - ncut*numpy.median(PS[:nband,:]) - c_pink
if c_pink<0:
  c_pink = 0.
else:
  c_pink = numpy.sqrt(c_pink)
if u_pink<0:
  u_pink = 0.
else:
  u_pink = numpy.sqrt(u_pink)

# reference output information
if refout:
  PS_plus  = numpy.sum(PS[nband+1,1:ncut+1]) - ncut*numpy.median(PS[nband+1,:])
  PS_minus = numpy.sum(PS[nband+2,1:ncut+1]) - ncut*numpy.median(PS[nband+2,:])
  # now need to get:
  # r_pink = amplitude of reference output 1/f noise
  # rho = correlation coefficient with common noise
  r_pink = numpy.sqrt(numpy.clip(PS_plus + PS_minus - 2*c_pink**2,0,None)/2.)
  rho = (PS_plus - PS_minus)/(4*r_pink*c_pink + 1e-24) # prevent division by zero
  rho = numpy.clip(rho,-1.,1.)
  #
  # and now the slope m (factor by which to multiply the correlated noise
  # to get what we see in the reference output)
  m_pink = r_pink*rho/(c_pink+1e-24)
  ru_pink = r_pink*(1-rho**2)

# normalize to variance (DN^2) per ln frequency
scale = numpy.sqrt(2/numpy.log(4096))
u_pink *= scale
c_pink *= scale
print('c_pink =', c_pink, 'u_pink =', u_pink)
if refout:
  ru_pink *= scale
  print('m_pink = ', m_pink, 'ru_pink =', ru_pink)

print('>>', c_pink**2, PS_plus, PS_minus)
for jj in range(numpy.shape(PS)[0]):
  print(numpy.sum(PS[jj,1:ncut+1]) - ncut*numpy.median(PS[jj,:]))

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
outcube = numpy.zeros((10,nside,nside))
outcube[0,:,:] = 65535 - medarray
outcube[1,:,:] = iqrarray/1.34896
outcube[2,:,:] = cdsnoise
outcube[3,:,:] = PCA0
outcube[4,:,:] = dark1
outcube[5,:,:] = dark2
outcube[6,:,:] = tnoise
# low CDS, high total noise pixels
outcube[7,:,:] = numpy.where(numpy.logical_and(cdsnoise<cds_cut,tnoise>numpy.median(tnoise)), 1, 0)
outcube[7,:4,:] = 0.; outcube[7,-4:,:] = 0.; outcube[7,:,:4] = 0.; outcube[7,:,-4:] = 0.
outcube[8,:,:] = dark1err
outcube[9,:,:] = dark2err
# label the slices
outcubeslices = ['BIAS', 'RESET', 'CDS', 'PCA0', 'DARK1', 'DARK2', 'TNOISE', 'LCDSHTN', 'DARK1ERR', 'DARK2ERR']

# array information on low CDS, high total noise pixels
badmap = outcube[7,4:-4,4:-4]
nbad1 = numpy.sum(badmap>.5)
nbad3 = numpy.sum(convolve(badmap,numpy.ones((3,3)),mode='same')>.5)
nbad5 = numpy.sum(convolve(badmap,numpy.ones((5,5)),mode='same')>.5)

primary_hdu = fits.PrimaryHDU()
out_hdu = fits.ImageHDU(outcube)
out_hdu.header['EXTNAME'] = 'NOISE'
out_hdu.header['TGROUP'] = tgroup
out_hdu.header['TDARK1'] = tgroup
out_hdu.header['TDARK2'] = tgroup*(ntslice-t_init)
out_hdu.header['ACN'] = acnnoise
out_hdu.header['C_PINK'] = c_pink
out_hdu.header['U_PINK'] = u_pink
out_hdu.header['PCA0_AMP'] = pca0_stdev
for k in range(nfile): out_hdu.header['NR_MF{:03d}'.format(k+1)] = filelist[k].strip()
# median information
out_hdu.header['CDS_CUT'] = (cds_cut, 'DN')
out_hdu.header['CDS_MED'] = (numpy.median(cdsnoise), 'DN')
out_hdu.header['TOT_MED'] = (numpy.median(tnoise), 'DN')
out_hdu.header['LCHTN1'] = (nbad1, 'low CDS, high total noise pixels')
out_hdu.header['LCHTN3'] = (nbad3, 'low CDS, high total noise pixels, 3x3 expanded')
out_hdu.header['LCHTN5'] = (nbad5, 'low CDS, high total noise pixels, 5x5 expanded')
out_hdu.header['NR_DATE'] = time.asctime(time.localtime(time.time()))

for idx in range(numpy.shape(outcube)[0]):
  out_hdu.header[outcubeslices[idx]] = (idx, 'Index of '+outcubeslices[idx])

ps_hdu = fits.ImageHDU(PS)

# CDS vs total noise 2D plot
d_np = .25
N_np = 80
ctnplot = numpy.zeros((N_np,N_np))
for j in range(N_np):
  tnmin = j*d_np; tnmax = tnmin + d_np
  if j==0: tnmin = -1e49
  if j==N_np-1: tnmax = 1e49
  for i in range(N_np):
    cnmin = i*d_np; cnmax = cnmin + d_np
    if i==0: cnmin = -1e49
    if i==N_np-1: cnmax = 1e49
    ctnplot[j,i] = numpy.count_nonzero(numpy.logical_and( numpy.logical_and(tnoise>=tnmin,tnoise<tnmax) ,
      numpy.logical_and(cdsnoise>=cnmin,cdsnoise<cnmax) ))
plot_hdu = fits.ImageHDU(ctnplot)
plot_hdu.header['INFO'] = '2D plot, CDS vs total noise'
plot_hdu.header['XDEF'] = ('CDS noise', 'DN')
plot_hdu.header['YDEF'] = ('total noise', 'DN')
plot_hdu.header['DNOISE'] = (d_np, 'DN')
plot_hdu.header['MAXNOISE'] = (d_np*N_np, 'DN')

outhdus = [primary_hdu, out_hdu, ps_hdu, plot_hdu]

# extra HDU if we asked for the reference output
if refout:

  extra = numpy.zeros((2,nside,nside//nch))

  xmin = nside
  xmax = nside + nside//nch
  nstep = ntslice - 2
  bandData = numpy.zeros((nfile2,nstep,nside,nside//nch), dtype=numpy.float32)
  for k in range(nfile2):
    for kb in range(nstep):
      bandData[k,kb,:,:] = pyirc.load_segment(filelist[k], fileformat, [xmin,xmax,0,nside], [t_init+kb  ], False).astype(numpy.float64)
  extra[0,:,:] = 65535-numpy.median(bandData, axis=(0,1))
  extra[1,:,:] = ( numpy.percentile(bandData, 75, axis=(0,1)) - numpy.percentile(bandData, 25, axis=(0,1)) ) / 1.34896

  amp33_hdu = fits.ImageHDU(extra.astype(numpy.float32))
  amp33_hdu.header['M_PINK'] = m_pink
  amp33_hdu.header['RU_PINK'] = ru_pink
  amp33_hdu.header['EXTNAME'] = 'AMP33'
  outhdus.append(amp33_hdu)

hdul = fits.HDUList(outhdus)
hdul.writeto(outfile, overwrite=True)

