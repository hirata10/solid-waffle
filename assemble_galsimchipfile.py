import sys
import time
import numpy
import re
from astropy.io import fits

# Basic parameters
noise_frame_bias = 0
noise_frame_ktc_plus_cds = 1
noise_frame_cds = 2
noise_frame_pca0 = 3

# Input file
if (len(sys.argv)<2):
  print('Needs configuration file.')
  exit()
config_file = sys.argv[1]
with open(config_file) as myf: content = myf.read().splitlines()
#
# Read information from configuration file
class EmptyClass:
  pass
configInfo = EmptyClass()
for line in content:
  #
  # Label
  m = re.search(r'^LABEL\:\s*(\S*)', line)
  if m: configInfo.label = m.group(1)
  m = re.search(r'^SCA\:\s*(\d+)', line)
  if m: configInfo.sca = int(m.group(1))
  m = re.search(r'^IN\:\s*(\S*)', line)
  if m: configInfo.IN = m.group(1)
  m = re.search(r'^OUT\:\s*(\S*)', line)
  if m: configInfo.OUT = m.group(1)
  m = re.search(r'^NOISE\:\s*(\S*)', line)
  if m: configInfo.NOISE = m.group(1)

# Check for information being there
if not hasattr(configInfo, 'NOISE'):
  print('Error: need NOISE')
  exit()

# Get information from summary file
# data
summaryData = numpy.loadtxt(configInfo.IN)
# summary information
with open(configInfo.IN) as myf: summaryinfo = myf.read().splitlines()
summaryMetadata = []
colData = []; countCol = False
for line in summaryinfo:
  m = re.search(r'^\#\ ', line)
  if m:
    mm = re.search(r'^\ *\d+\,(.*)', line[2:])
    if mm: colData.append(mm.group(1))
    if line[2:] == r'Columns\:': countCol = True
    if not countCol: summaryMetadata.append(line[2:])

nx = 1 + int(numpy.max(summaryData[:,0]))
ny = 1 + int(numpy.max(summaryData[:,1]))
print('Main file grid: nx = {:d}, ny = {:d}'.format(nx,ny))

# General map
gaincol=-1
for k in range(len(colData)):
  m = re.search(r'alpha\,beta\-corrected gain', colData[k])
  if m: gaincol=k
if gaincol>=0:
  print('Gain column -> {:d}'.format(gaincol))
else:
  print('Gain column not found.')
  exit()
gain = summaryData[:,gaincol].reshape((ny,nx))
#
# ... and here, SuperPixelMask tells us which super-pixels were masked
SuperPixelMask = numpy.where(gain<1e-99)
nmask = numpy.size(SuperPixelMask[0])
print('SuperPixelMask:', SuperPixelMask)
print('.. length:', nmask)
#
# ... we're going to replace masked data with the median of the array
medgain = numpy.median(gain)
gain[SuperPixelMask] = medgain

# Output information
primary_hdu = fits.PrimaryHDU()
primary_hdu.header['LABEL'] = configInfo.label
primary_hdu.header['SCA'] = configInfo.sca
primary_hdu.header['INFILE'] = (configInfo.IN.strip(), 'Main input')
primary_hdu.header['IN_NOISE'] = (configInfo.NOISE.strip(), 'Noise input')
primary_hdu.header['OUTFILE'] = configInfo.OUT.strip()
primary_hdu.header['GENDATE'] = time.asctime(time.localtime(time.time()))
primary_hdu.header['SR_NY'] = (ny, 'ny for summary file')
primary_hdu.header['SR_NX'] = (nx, 'nx for summary file')
primary_hdu.header['SR_NMASK'] = (nmask, 'number of masked super pixels')

# Source information
src_hdu = fits.ImageHDU()
src_hdu.header['EXTNAME'] = 'SOURCES'
#
# copy metadata
for k in range(len(summaryMetadata)):
  src_hdu.header['META{:d}'.format(k)] = summaryMetadata[k].strip()
for k in range(len(colData)):
  src_hdu.header['INCOL{:d}'.format(k)] = colData[k].strip()

# Read noise HDU
noisefile = fits.open(configInfo.NOISE)
#
# CDS noise data
cdsnoise = noisefile['NOISE'].data[noise_frame_cds,:,:]
med_cdsnoise_DN = numpy.median(cdsnoise[4:-4,4:-4])
med_cdsnoise_DN_ref = numpy.median(numpy.concatenate((cdsnoise[:,:4].flatten(), cdsnoise[:,-4:].flatten(),\
  cdsnoise[:4,4:-4].flatten(), cdsnoise[-4:,4:-4].flatten())))
print('CDS noise (DN): median ref, median active', med_cdsnoise_DN_ref, med_cdsnoise_DN)
if med_cdsnoise_DN_ref>med_cdsnoise_DN:
  print('Warning: med_cdsnoise_DN_ref>med_cdsnoise_DN')
# reset noise
med_ktccdsnoise_DN = numpy.median(noisefile['NOISE'].data[noise_frame_ktc_plus_cds,4:-4,4:-4])
if med_ktccdsnoise_DN>med_cdsnoise_DN:
  med_ktcnoise_DN = numpy.sqrt(med_ktccdsnoise_DN**2-med_cdsnoise_DN**2)
  print('Reset noise (DN):', med_ktcnoise_DN)
else:
  med_ktcnoise_DN = 0.
  print('Warning: med_ktccdsnoise_DN<=med_cdsnoise_DN')
#
for noisekey in noisefile['NOISE'].header.keys():
  m = re.search(r'^NR\_MF', noisekey)
  if m:
    src_hdu.header[noisekey] = (noisefile['NOISE'].header[noisekey], 'noise file')
noise_hdu = fits.ImageHDU(noisefile['NOISE'].data[noise_frame_pca0,:,:]) # pca0 in last slice
noise_hdu.header['EXTNAME'] = 'READ'
# additional keywords
noise_hdu.header['NGNAXIS1'] = 4096
noise_hdu.header['NGNAXIS2'] = 4096
noise_hdu.header['NG_NOUT'] = 32
#
# these are default values, we may update them later
# also they don't have guide window information
noise_hdu.header['NG_NFOH'] = (1, 'placeholder')
noise_hdu.header['NG_NROH'] = (12, 'placeholder')
#
noise_hdu.header['DT'] = (5e-6, '200 kHz sampling')
#
# more noise parameters, placeholders where indicated
# note noise properties are in e
noise_hdu.header['PEDESTAL'] = (4., 'placeholder')
noise_hdu.header['RD_NOISE'] = med_cdsnoise_DN*numpy.median(gain) # convert from DN --> e
noise_hdu.header['C_PINK'] = (3., 'placeholder')
noise_hdu.header['U_PINK'] = (1., 'placeholder')
noise_hdu.header['ACN'] = (0.5, 'placeholder')
noise_hdu.header['PCA0_AMP'] = noisefile['NOISE'].header['PCA0_AMP']*numpy.median(gain) # convert from DN --> e
noise_hdu.header['REFPIXNR'] = med_cdsnoise_DN_ref/med_cdsnoise_DN
noise_hdu.header['KTCNOISE'] = med_ktcnoise_DN*numpy.median(gain) # convert from DN --> e

noise_hdu.header.add_comment('Noise properties in electrons, not DN')

# Make gain map
gain_hdu= fits.ImageHDU(gain)
gain_hdu.header['EXTNAME'] = 'GAIN'

# Final output step
hdul = fits.HDUList([primary_hdu, src_hdu, noise_hdu, gain_hdu])
hdul.writeto(configInfo.OUT, overwrite=True)
