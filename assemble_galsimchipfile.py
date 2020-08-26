import sys
import time
import numpy
import re
import pyirc
from astropy.io import fits
from scipy import ndimage

#######################################################################
# Setup
#######################################################################

# Basic parameters
noise_frame_bias = 0
noise_frame_ktc_plus_cds = 1
noise_frame_cds = 2
noise_frame_pca0 = 3
noise_frame_dark1 = 4
noise_frame_dark2 = 5

# Choice of input format
informat = 4
d_offset = 1 # 1 if first frame is *after* the reset, 0 if it is reset-read

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
configInfo.use_spr = False
configInfo.flatBFE = False
configInfo.FW = False
configInfo.Dark = False
configInfo.PersistScript = ''
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
  m = re.search(r'^SPR\:\s*(\S*)\+(\d+),(\d+)', line)
  if m:
    configInfo.use_spr = True
    configInfo.SPR = m.group(1)
    configInfo.NSPR = int(m.group(2))
    configInfo.dSPR = int(m.group(3))
  m = re.search(r'^FLATBFE', line)
  if m: configInfo.flatBFE = True
  m = re.search(r'^FW\:\s*(\S*)\+(\d+)\s+(\d+)$', line)
  if m:
    configInfo.NLD = int(m.group(3)) # non-linear order D
    configInfo.NFullWellFile = int(m.group(2))
    configInfo.FullWellFile = m.group(1)
    configInfo.FW = True
  m = re.search(r'^DARK\:\s*(\S*)\+(\d+)$', line)
  if m:
    configInfo.NDarkFile = int(m.group(2))
    configInfo.DarkFile = m.group(1)
    configInfo.Dark = True
  m = re.search(r'^PERSIST\:\s*(\S*)', line)
  if m: configInfo.PersistScript = m.group(1)

# Check for information being there
if not hasattr(configInfo, 'NOISE'):
  print('Error: need NOISE')
  exit()

#######################################################################
# Get information from the summary file
#######################################################################

badpix = numpy.zeros((4096,4096), dtype=numpy.uint16)

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

col = pyirc.IndexDictionary(0)
sbfe = 2
col.addbfe(sbfe)
p = len(colData)-col.N-2; print('number of nl coefficients', p)
col.addhnl(p)

# General map
gain = summaryData[:,col.g+2].reshape((ny,nx))
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
primary_hdu.header['ELECTR'] = ('Lab', 'bias/gain/cnl/noise from lab electronics, flight will be different')

#
# copy metadata
this_metadata = []
for k in range(len(summaryMetadata)):
  this_metadata.append(summaryMetadata[k])
for k in range(len(colData)):
  this_metadata.append(colData[k])

# Get BFE information
ssbfe = 2*sbfe+1
bfeData = numpy.zeros((ssbfe,ssbfe,ny,nx))
for x in range(ssbfe):
  for y in range(ssbfe):
    thisdata = numpy.copy(summaryData[:,col.Nb+ssbfe*y+x+2].reshape((ny,nx)))
    if configInfo.flatBFE:
      thisdata[:,:] = numpy.median(thisdata)
    else:
      thisdata[SuperPixelMask] = numpy.median(thisdata)
    bfeData[y,x,:,:] = thisdata
bfe_hdu = fits.ImageHDU(bfeData.astype(numpy.float32))
bfe_hdu.header['EXTNAME'] = 'BFE'

#######################################################################
# IPC & VTPE HDUs
#######################################################################

# Get IPC information
ipc_grid = 8
nside = 4096
nipc = nside//ipc_grid
ipc_auto = numpy.zeros((3,3,nipc,nipc))
#
# auto-correlation map
aH = summaryData[:,col.alphaH+2].reshape((ny,nx))
aV = summaryData[:,col.alphaV+2].reshape((ny,nx))
aD = summaryData[:,col.alphaD+2].reshape((ny,nx))
aH[SuperPixelMask] = numpy.median(aH)
aV[SuperPixelMask] = numpy.median(aV)
aD[SuperPixelMask] = numpy.median(aD)
dy = nipc//ny; dx = nipc//nx
for j in range(ny):
  for i in range(nx):
    ipc_auto[1,2,j*dy:j*dy+dy,i*dx:i*dx+dx] = aH[j,i]
    ipc_auto[2,1,j*dy:j*dy+dy,i*dx:i*dx+dx] = aV[j,i]
    ipc_auto[2,2,j*dy:j*dy+dy,i*dx:i*dx+dx] = aD[j,i]
# copy by symmetry
ipc_auto[1,0,:,:] = ipc_auto[1,2,:,:]
ipc_auto[0,1,:,:] = ipc_auto[2,1,:,:]
ipc_auto[2,0,:,:] = ipc_auto[2,2,:,:]
ipc_auto[0,2,:,:] = ipc_auto[2,2,:,:]
ipc_auto[0,0,:,:] = ipc_auto[2,2,:,:]
# normalize
ipc_auto[1,1,:,:] = 0.
sum_ipc = numpy.zeros((nipc,nipc))
for dj in range(3):
  for di in range(3):
    sum_ipc += ipc_auto[dj,di,:,:]
ipc_auto[1,1,:,:] =  1.-sum_ipc
ipc_full = numpy.copy(ipc_auto)

# Placeholder for VTPE
vtpe_good = False
vtpe = numpy.zeros((3,512,512)); vtpe[2,:,:] = 1.

# update with SPR data, if available
if configInfo.use_spr:
  SPRFile = []; medsig = []
  ipc_full = numpy.zeros((3,3,512,512))
  sx = 512//nx; sy = 512//ny
  for iy in range(ny):
    for ix in range(nx):
      for ky in range(3):
        for kx in range(3):
          ipc_full[ky,kx,iy*sy:(iy+1)*sy,ix*sx:(ix+1)*sx] = ipc_auto[ky,kx,iy,ix]
  alphadata = numpy.zeros((configInfo.NSPR, 13, 512, 512))
  m = re.search(r'^(.+)(\d{2})\_alpha\.fits', configInfo.SPR)
  if not m:
    print('Error: pattern match failed on ', configInfo.SPR)
    exit()
  st = m.group(1); en = int(m.group(2))
  for j in range(configInfo.NSPR):
    SPRFile.append(st + '{:02d}_alpha.fits'.format(en))
    this_metadata.append('Used SPR file: '+SPRFile[j])
    with fits.open(SPRFile[j]) as G:
      for keyword in G[0].header.keys():
        m = re.search(r'^ARGV', keyword)
        if m: this_metadata.append('  '+keyword+' -> '+G[0].header[keyword])
        m = re.search(r'^INF', keyword)
        if m: this_metadata.append('  '+keyword+' -> '+G[0].header[keyword])
      medsig.append(float(G[0].header['MEDSIG']))
      alphadata[j,:,:,:] = G[0].data
    print(' median signal:', medsig[j], 'DN'); this_metadata.append(' median signal: '+str(medsig[j])+' DN')
    en+=configInfo.dSPR
  medsig = numpy.asarray(medsig)

  # put into IPC maps
  this_aV = numpy.median(alphadata[-3:,9,:,:], axis=0)
  this_aH = numpy.median(alphadata[-3:,0,:,:], axis=0)
  this_aLL = numpy.median(alphadata[-3:,8,:,:], axis=0)
  this_aLR = numpy.median(alphadata[-3:,10,:,:], axis=0)
  s = 9
  s_th = .01 # threshold for clipping
  this_aV2 = ndimage.median_filter(this_aV, size=s, mode='reflect')
  this_aV = numpy.where(numpy.absolute(this_aV2-this_aV)>s_th, this_aV2, this_aV)
  this_aH2 = ndimage.median_filter(this_aH, size=s, mode='reflect')
  this_aH = numpy.where(numpy.absolute(this_aH2-this_aH)>s_th, this_aH2, this_aH)
  s_th = .001 # threshold for clipping
  this_aLL2 = ndimage.median_filter(this_aLL, size=s, mode='reflect')
  this_aLL = numpy.where(numpy.absolute(this_aLL2-this_aLL)>s_th, this_aLL2, this_aLL)
  this_aLR2 = ndimage.median_filter(this_aLR, size=s, mode='reflect')
  this_aLR = numpy.where(numpy.absolute(this_aLR2-this_aLR)>s_th, this_aLR2, this_aLR)
  del this_aV2; del this_aH2; del this_aLL2; del this_aLR2

  ipc_full[2,1,:,:] = ipc_full[0,1,:,:] = this_aV
  ipc_full[1,2,:,:] = ipc_full[1,0,:,:] = this_aH
  ipc_full[2,2,:,:] = ipc_full[0,0,:,:] = this_aLL
  ipc_full[2,0,:,:] = ipc_full[0,2,:,:] = this_aLR

  # central pixel correction
  sum_ipc = numpy.zeros((512,512))
  ipc_full[1,1,:,:] = 0.
  for dj in range(3):
    for di in range(3):
      sum_ipc += ipc_full[dj,di,:,:]
  ipc_full[1,1,:,:] =  1.-sum_ipc

  # VTPE maps
  vt_big = alphadata[-1,9,:,:] - alphadata[-1,5,:,:]
  vt_small = alphadata[0,9,:,:] - alphadata[0,5,:,:]
  #vt_big = numpy.median(alphadata[-3:,9,:,:] - alphadata[-3:,5,:,:], axis=0)
  #vt_small = numpy.median(alphadata[:3,9,:,:] - alphadata[:3,5,:,:], axis=0)
  s_th = .005
  vt_big2 = ndimage.median_filter(vt_big, size=s, mode='reflect')
  vt_big = numpy.where(numpy.absolute(vt_big-vt_big2)>s_th, vt_big2, vt_big)
  vt_small2 = ndimage.median_filter(vt_small, size=s, mode='reflect')
  vt_small = numpy.where(numpy.absolute(vt_small-vt_small2)>s_th, vt_small2, vt_small)
  g_copy = numpy.zeros((512,512))
  for iy in range(ny):
    for ix in range(nx):
      g_copy[iy*sy:(iy+1)*sy,ix*sx:(ix+1)*sx] = gain[iy,ix]
  Q_big = medsig[-1] * g_copy
  Q_small = medsig[0] * g_copy
  vtpe[1,:,:] = numpy.minimum((vt_big-vt_small)/numpy.log(Q_big/Q_small), -1.e-24*numpy.ones((512,512), dtype=numpy.float64)) # to guarantee not zero
  vtpe[1,:,:] = ndimage.median_filter(vtpe[1,:,:], size=(5,1), mode='reflect')
  #
  # floor
  dy = nipc//ny; dx = nipc//nx
  thisfloor = numpy.zeros((ny,nx))
  for iy in range(ny):
    for ix in range(nx):
      thisfloor[iy,ix] = 2*(numpy.average(this_aV[iy*dy:(iy+1)*dy,ix*dx:(ix+1)*dx])-aV[iy,ix])
  thisfloor = ndimage.median_filter(thisfloor, size=(3,1), mode='reflect')
  for iy in range(ny):
    for ix in range(nx):
      xicpt = numpy.log(numpy.average(Q_small[iy*dy:(iy+1)*dy,ix*dx:(ix+1)*dx]))\
        - numpy.average(vt_small[iy*dy:(iy+1)*dy,ix*dx:(ix+1)*dx]-thisfloor[iy,ix])/numpy.average(vtpe[1,iy*dy:(iy+1)*dy,ix*dx:(ix+1)*dx])
      print(iy,ix,thisfloor[iy,ix],numpy.average(Q_small[iy*dy:(iy+1)*dy,ix*dx:(ix+1)*dx]),numpy.average(vt_small[iy*dy:(iy+1)*dy,ix*dx:(ix+1)*dx]),\
        numpy.average(vtpe[1,iy*dy:(iy+1)*dy,ix*dx:(ix+1)*dx]),xicpt)
      vtpe[2,iy*dy:(iy+1)*dy,ix*dx:(ix+1)*dx] = 1.+numpy.exp(xicpt)
  vtpe[0,:,:] = vt_big - vtpe[1,:,:]*numpy.log(1. + Q_big/vtpe[2,:,:])

  del alphadata # cleanup

del ipc_auto
del sum_ipc

# write IPC information
ipc_hdu = fits.ImageHDU(ipc_full.astype(numpy.float32))
ipc_hdu.header['EXTNAME'] = 'IPC'

(ay,ax,nyi,nxi)=numpy.shape(ipc_full)
ipcflat_hdu = fits.ImageHDU(ipc_full.astype(numpy.float32).reshape((ay*ax,nyi,nxi)))
ipcflat_hdu.header['EXTNAME'] = 'IPCFLAT'

# VTPE information
vtpe_hdu = fits.ImageHDU(vtpe.astype(numpy.float32))
vtpe_hdu.header['EXTNAME'] = 'VTPE'
vtpe_hdu.header['SLICE01'] = ('MAP OF A_VTPE', 'floor, dimensionless')
vtpe_hdu.header['SLICE02'] = ('MAP OF B_VTPE', 'slope, dimensionless')
vtpe_hdu.header['SLICE03'] = ('MAP OF DQ0_VTPE', 'break in electrons')
if vtpe_good:
  vtpe_hdu.header['ISGOOD'] = (True, 'Full VTPE information')
else:
  vtpe_hdu.header['ISGOOD'] = (False, '** VTPE IS PLACEHOLDER **')

#######################################################################
# Read noise HDU
#######################################################################
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
  this_metadata.append('Noise header:' + noisekey + '        ' + str(noisefile['NOISE'].header[noisekey]))
noise_hdu = fits.ImageHDU(noisefile['NOISE'].data[noise_frame_pca0,:,:])
noise_hdu.header['EXTNAME'] = 'READ'
# additional keywords
noise_hdu.header['NGNAXIS1'] = 4096
noise_hdu.header['NGNAXIS2'] = 4096
noise_hdu.header['NG_NOUT'] = 32
#
# these are default values, we may update them later
# also they don't have guide window information
noise_hdu.header['NG_NFOH'] = (1, 'placeholder')
noise_hdu.header['NG_NROH'] = (7, 'in lab')
#
noise_hdu.header['DT'] = (5e-6, '200 kHz sampling')
#
# more noise parameters, placeholders where indicated
# note noise properties are in e
noise_hdu.header['PEDESTAL'] = (4., 'placeholder')
noise_hdu.header['RD_NOISE'] = med_cdsnoise_DN*numpy.median(gain) # convert from DN --> e
noise_hdu.header['C_PINK'] = noisefile['NOISE'].header['C_PINK']*numpy.median(gain) # convert from DN --> e
noise_hdu.header['U_PINK'] = noisefile['NOISE'].header['U_PINK']*numpy.median(gain) # convert from DN --> e
noise_hdu.header['ACN'] = noisefile['NOISE'].header['ACN']*numpy.median(gain) # convert from DN --> e
noise_hdu.header['PCA0_AMP'] = noisefile['NOISE'].header['PCA0_AMP']*numpy.median(gain) # convert from DN --> e
noise_hdu.header['REFPIXNR'] = med_cdsnoise_DN_ref/med_cdsnoise_DN
noise_hdu.header['KTCNOISE'] = med_ktcnoise_DN*numpy.median(gain) # convert from DN --> e

noise_hdu.header.add_comment('Noise properties in electrons, not DN')

# Get the bias information
bias_hdu = fits.ImageHDU(numpy.clip(noisefile['NOISE'].data[noise_frame_bias,:,:], 0, 65535).astype(numpy.uint16))
bias_hdu.header['EXTNAME'] = 'BIAS'

#######################################################################
# Dark current HDU (uses some information from the noise for hot pixels)
#######################################################################

# dark current
td1 = noisefile['NOISE'].header['TDARK1']
td2 = noisefile['NOISE'].header['TDARK2']
dark_current = numpy.copy(noisefile['NOISE'].data[noise_frame_dark1,:,:])
dark_crit1 = 1024./td2; dark_crit2 = 0.
this_metadata.append('dark threshold 1: '+str(dark_crit1)+' DN/s')
badpix = numpy.bitwise_or(numpy.where(dark_current<dark_crit1, 0, 4).astype(badpix.dtype), badpix)
dark_current = numpy.where(dark_current<dark_crit1, noisefile['NOISE'].data[noise_frame_dark2,:,:], dark_current)
dark_var_rate = numpy.zeros((4096,4096), dtype=numpy.float32)

# if there are separate dark files, use them!
if configInfo.Dark:
  # first, get the files
  this_metadata.append('')
  this_metadata.append('Dark current computation:')
  DarkFile = []
  m = re.search(r'^(.+)(\d{3})\.fits', configInfo.DarkFile)
  if not m:
    print('Error: pattern match failed on ', configInfo.DarkFile)
    exit()
  st = m.group(1); en = int(m.group(2))
  for j in range(configInfo.NDarkFile):
    DarkFile.append(st + '{:03d}.fits'.format(en))
    this_metadata.append('Used dark file: '+DarkFile[j])
    en+=1
  nt = pyirc.get_num_slices(informat, DarkFile[0])
  NFowler = 4
  with fits.open(DarkFile[0]) as G:
    td3 = float(G[0].header['TGROUP'])*(nt-NFowler)
  darkimage = numpy.zeros((configInfo.NDarkFile, 4096, 4096), dtype=numpy.float32)   # Fowler 4
  for j in range(configInfo.NDarkFile):
    for k in range(NFowler):
      thisdiff = pyirc.load_segment(DarkFile[j], informat, [0,4096,0,4096], [1+k], False).astype(numpy.float32)\
                         - pyirc.load_segment(DarkFile[j], informat, [0,4096,0,4096], [nt-k], False).astype(numpy.float32)
      darkimage[j,:,:] = darkimage[j,:,:] + thisdiff
    darkimage[j,:,:] /= NFowler
    # reference pixel subtraction
    for i in range(4096):
      dside = numpy.median(numpy.concatenate((darkimage[j,i,:4],darkimage[j,i,-4:])))
      darkimage[j,i,:] = darkimage[j,i,:] - dside
    for i in range(32):
      xmin = 128*i; xmax = xmin+128
      dside = numpy.median(darkimage[j,-4:,xmin:xmax])
      darkimage[j,:,xmin:xmax] = darkimage[j,:,xmin:xmax] - dside
  longdark = numpy.median(darkimage, axis=0).astype(numpy.float32)/td3
  dark_var_rate = (numpy.median(numpy.absolute(darkimage/td3-longdark), axis=0).astype(numpy.float32)/0.67448)**2*td3
  # substitute this dark rate if need be
  dark_crit2 = 256./td2
  this_metadata.append('dark threshold 2: '+str(dark_crit2)+' DN/s')
  badpix = numpy.bitwise_or(numpy.where(dark_current<dark_crit2, 0, 2).astype(badpix.dtype), badpix)
  dark_current = numpy.where(dark_current<dark_crit2, longdark, dark_current)
  del darkimage

# gain conversion: convert dark from DN/s -> e/s
for j in range(ny):
  ymin = j*4096//ny; ymax = ymin+4096//ny
  for i in range(nx):
    xmin = i*4096//nx; xmax = xmin+4096//nx
    dark_current[ymin:ymax,xmin:xmax] *= gain[j,i]
    dark_var_rate[ymin:ymax,xmin:xmax] *= gain[j,i]**2

dark_hdu = fits.ImageHDU(dark_current.astype(numpy.float32))
dark_hdu.header['EXTNAME'] = 'DARK'
dark_hdu.header['DARK_CR1'] = (dark_crit1, 'DN/s')
dark_hdu.header['DARK_CR2'] = (dark_crit2, 'DN/s')
dark_hdu.header['TDARK1'] = (td1, '1st tier dark time in seconds')
dark_hdu.header['TDARK2'] = (td2, '2nd tier dark time in seconds')
dark_hdu.header['TDARK3'] = (td3, '3rd tier dark time in seconds')
dark_hdu.header['LONGDARK'] = (configInfo.Dark, 'Long darks available?')
dark_hdu.header['LDFOWLER'] = (NFowler, 'Fowler-n sampling of long dark')

darkvar_hdu = fits.ImageHDU(dark_var_rate.astype(numpy.float32))
darkvar_hdu.header['EXTNAME'] = 'DARKVAR'
darkvar_hdu.header['LONGDARK'] = (configInfo.Dark, 'Long darks available?')
darkvar_hdu.header['LDFOWLER'] = (NFowler, 'Fowler-n sampling of long dark')


#######################################################################
# Non-linearity and saturation information
#######################################################################

# Extract full well information
if configInfo.FW:
  # first, get the files
  this_metadata.append('')
  this_metadata.append('Full well computation:')
  FWFile = []
  m = re.search(r'^(.+)(\d{3})\.fits', configInfo.FullWellFile)
  if not m:
    print('Error: pattern match failed on ', configInfo.FullWellFile)
    exit()
  st = m.group(1); en = int(m.group(2))
  for j in range(configInfo.NFullWellFile):
    FWFile.append(st + '{:03d}.fits'.format(en))
    this_metadata.append('Used full well file: '+FWFile[j])
    en+=1
  nt = pyirc.get_num_slices(informat, FWFile[0])
  print('Files:', FWFile[0], ' ... ', FWFile[-1], ', nt=', nt)
  with fits.open(FWFile[0]) as G: tgfw = float(G[0].header['TGROUP']) # get group time
  print('tgfw =', tgfw, 's')
  my_stack = numpy.zeros((nt, 4096, 4096))
  tempstack = numpy.zeros((configInfo.NFullWellFile, 4096, 4096), dtype=numpy.float32)
  for i in range(nt):
    for j in range(configInfo.NFullWellFile):
      tempstack[j,:,:] = pyirc.load_segment(FWFile[j], informat, [0,4096,0,4096], [1], False).astype(numpy.float32)\
                         - pyirc.load_segment(FWFile[j], informat, [0,4096,0,4096], [i+1], False).astype(numpy.float32)
    my_stack[i,:,:] = numpy.median(tempstack, axis=0).astype(numpy.float64)
      # convert to double precision for polynomial fitting
    print('time step {:2d} --> median signal (rel. to first) = {:9.2f}'.format(i+1, numpy.median(my_stack[i,:,:])),\
      time.asctime(time.localtime(time.time()))); sys.stdout.flush()
  del tempstack
  
  # allocate arrays for poly coefficients, full well, time stamps
  tmax = numpy.zeros((4096,4096), dtype=numpy.int16); tmax[:,:] = nt-1
  poly_coefs = numpy.zeros((configInfo.NLD+1, 4096, 4096))
  these_poly_coefs = numpy.zeros((configInfo.NLD+1, 4096, 4096))
  #
  # now do fit over each range
  print('using order', configInfo.NLD)
  for tm in range(configInfo.NLD+1,nt)[::-1]:
    A = numpy.zeros((configInfo.NLD+1,configInfo.NLD+1))
    B = numpy.zeros((configInfo.NLD+1,tm))
    for i in range(configInfo.NLD+1):
      for j in range(configInfo.NLD+1):
        A[i,j] = numpy.sum(numpy.array(range(d_offset,tm+d_offset)).astype(numpy.float64)**(i+j))
      for j in range(tm):
        B[i,j] = (d_offset+j)**i
    AinvB = numpy.matmul(numpy.linalg.inv(A), B)
    #print('shapes: ', numpy.shape(AinvB), numpy.shape(my_stack))
    these_poly_coefs[:,:,:] = numpy.tensordot(AinvB, my_stack[:tm,:,:], axes=([1],[0])) # sum_j Ainv[i,j] * my_stack[j,:,:]
    #
    if tm==nt-1: poly_coefs[:,:,:] = these_poly_coefs # start by accepting all coefficients
    #
    # if there is something wrong with the previous fit, accept this one
    if tm<nt-1:
      err = numpy.copy(my_stack[:tm,:,:])
      t = numpy.asarray(range(d_offset, d_offset+tm)).astype(numpy.float64)
      for i in range(configInfo.NLD+1): err -= numpy.tensordot(t**i, these_poly_coefs[i,:,:], axes=0)
      err = numpy.sqrt(numpy.average(err**2, axis=0))
      BadMask = numpy.logical_or( my_stack[tm,:,:]-my_stack[tm-1,:,:]<(my_stack[tm-1,:,:]-my_stack[tm-2,:,:])/2., err>327.68)
      der = numpy.zeros((4096,4096))
      for i in range(1,configInfo.NLD+1): der += poly_coefs[i,:,:]*(d_offset+tmax-1).astype(numpy.float64)**(i-1)*i
      BadMask = numpy.logical_or(BadMask, der<0)
      print('tm = ', tm, ' array-median RMS err = ', numpy.median(err), 'mask=', numpy.sum(numpy.where(BadMask,1,0))); sys.stdout.flush()
      poly_coefs = numpy.where(BadMask, these_poly_coefs, poly_coefs)
      tmax = numpy.where(BadMask, tm, tmax)
  del these_poly_coefs
  #
  for tm in range(1,configInfo.NLD+1)[::-1]:
    der = numpy.zeros((4096,4096))
    for i in range(1,configInfo.NLD+1): der += poly_coefs[i,:,:]*float(d_offset+tm)**(i-1)*i
    delta = numpy.zeros((4096,4096))
    for i in range(1,configInfo.NLD+1): delta += poly_coefs[i,:,:]*(float(d_offset+tm)**i - float(d_offset+tm-1)**i)
    tmax = numpy.where(numpy.logical_or(der<0,delta<0), tm, tmax)
  del der; del delta
  # print statistics of the tmax file
  for tm in range(0,nt): print('{:2d} {:8d}'.format(tm, numpy.sum(numpy.where(tmax==tm,1,0).astype(numpy.int32))))
  # figure out saturation level
  sat_level = poly_coefs[1,:,:]*(d_offset+tmax-1)
  dx = 4096//nx; dy = 4096//ny
  for iy in range(ny):
    for ix in range(nx):
      sat_level[dy*iy:dy*(iy+1),dx*ix:dx*(ix+1)] *= gain[iy,ix] # convert to electrons
  sat_level[:4,:] = 0.; sat_level[-4:,:] = 0.; sat_level[:,:4] = 0.; sat_level[:,-4:] = 0.
  print('saturation', sat_level[4:-4:113,4:-4:113]) # sample every 113th in x,y for display
  # make error map
  err_level = numpy.zeros((4096,4096), dtype=numpy.float32)
  for tm in range(1,nt):
    err = numpy.copy(my_stack[:tm,:,:])
    t = numpy.asarray(range(d_offset, d_offset+tm)).astype(numpy.float64)
    for i in range(configInfo.NLD+1): err -= numpy.tensordot(t**i, poly_coefs[i,:,:], axes=0)
    err_level = numpy.where(tmax==tm, numpy.sqrt(numpy.mean(err**2,axis=0)).astype(numpy.float32), err_level)
  #
  # convert to output cube
  # (and zero out reference pixels)
  poly_coefs[:,:4,:] = 0.; poly_coefs[:,-4:,:] = 0.; poly_coefs[:,:,:4] = 0.; poly_coefs[:,:,-4:] = 0.
  poly_coefs[1,:4,:] = 1.; poly_coefs[1,-4:,:] = 1.; poly_coefs[1,:,:4] = 1.; poly_coefs[1,:,-4:] = 1.
  for q in range(1,configInfo.NLD-1):
    poly_coefs[1+q,:,:] = -poly_coefs[1+q,:,:]/poly_coefs[1,:,:]**(1+q) # write in terms of DN_lin
    dx = 4096//nx; dy = 4096//ny
    for iy in range(ny):
      for ix in range(nx):
        poly_coefs[1+q,dy*iy:dy*(iy+1),dx*ix:dx*(ix+1)] /= gain[iy,ix]**q # convert to electrons
  flat_field = numpy.copy(poly_coefs[1,:,:])
  for iy in range(ny):
    for ix in range(nx):
      flat_field[dy*iy:dy*(iy+1),dx*ix:dx*(ix+1)] *= gain[iy,ix]
  flat_field[:4,:] = 0.; flat_field[-4:,:] = 0.; flat_field[:,:4] = 0.; flat_field[:,-4:] = 0. # take out reference pixels
  flat_field[4:-4,4:-4] -= dark_current[4:-4,4:-4]*tgfw # take out the dark
  flat_field /= numpy.median(flat_field)
  flat_field = numpy.maximum(flat_field,0)
  poly_coefs = poly_coefs[2:,:,:]
# alternative
else:
  this_metadata.append('')
  this_metadata.append('Full well computations: ***PLACEHOLDER***')
  sat_level = numpy.zeros((4096,4096)) + 80000 # placeholder
  poly_coefs = numpy.zeros((p,ny,nx))
  for q in range(p): poly_coefs[q,:,:] = numpy.copy(summaryData[:,col.Nbb+2:col.Nbb+q+2].reshape((ny,nx)))/gain**q
  flat_field = numpy.copy(summaryData[:,col.I].reshape((ny,nx)))
  err_level = numpy.zeros((4096,4096), dtype=numpy.float32)
#
# make CNL + saturation HDUs
cnl_hdu = fits.ImageHDU(poly_coefs)
cnl_hdu.header['EXTNAME'] = 'CNL'
cnl_hdu.header['ERR50'] = (numpy.percentile(err_level, 50.), '50th percentile CNL fit error (DN)')
cnl_hdu.header['ERR90'] = (numpy.percentile(err_level, 90.), '90th percentile CNL fit error (DN)')
cnl_hdu.header['ERR95'] = (numpy.percentile(err_level, 95.), '95th percentile CNL fit error (DN)')
cnl_hdu.header['ERR99'] = (numpy.percentile(err_level, 99.), '99th percentile CNL fit error (DN)')
saturate_hdu= fits.ImageHDU(numpy.floor(sat_level).astype(numpy.int32))
saturate_hdu.header['EXTNAME'] = 'SATURATE'
#
# relative QE
relqe1_hdu = fits.ImageHDU(flat_field.astype(numpy.float32))
relqe1_hdu.header['EXTNAME'] = 'RELQE1'
#
this_metadata.append('')

# Flag bad pixels from this test
badpix = numpy.bitwise_or(badpix, numpy.where(flat_field<0.25,1,0).astype(badpix.dtype))

# find weird pixels
flat2 = numpy.copy(flat_field)
dy = 4096//ny; dx = 4096//nx
for iy in range(ny):
  for ix in range(nx):
    flat2[dy*iy:dy*(iy+1),dx*ix:dx*(ix+1)] /= numpy.median(flat_field[dy*iy:dy*(iy+1),dx*ix:dx*(ix+1)])
weirdmap = numpy.where(numpy.absolute(flat2-1)>.1, 1, 0)
weirdmap[4:-4,4:-4] = ndimage.maximum_filter(weirdmap[4:-4,4:-4], size=3, mode='reflect')
badpix[4:-4,4:-4] = numpy.bitwise_or(badpix[4:-4,4:-4], numpy.where(weirdmap[4:-4,4:-4],8,0).astype(badpix.dtype))

#######################################################################
# Persistence information
#######################################################################

if configInfo.PersistScript:
  this_metadata.append('Persistence calculation:')
  with open(configInfo.PersistScript) as myf: content = myf.read().splitlines()
  NP = len(content)
  Q = []
  persistence_map = numpy.zeros((NP,4096,4096), dtype=numpy.float32)
  for f in range(NP):
    m = re.search(r'^([\d\.\+\-Ee]+)\s*(\S*)', content[f])
    if m:
      Q.append(float(m.group(1))); fn = m.group(2)
      this_metadata.append('Reading ... '+content[f])
      with fits.open(DarkFile[0]) as G:
        frtime = float(G[0].header['FRTIME'])
        tgroup = float(G[0].header['TGROUP'])
        ngroups = float(G[0].header['NGROUPS'])
      ngroups = min(4,ngroups) # stop after 4th group
      print('<--', fn); sys.stdout.flush()
      persistence_map[f,:,:] = pyirc.load_segment(fn, informat, [0,4096,0,4096], [1], False).astype(numpy.float32)
      for k in range(2,ngroups+1):
        persistence_map[f,:,:] = persistence_map[f,:,:] - pyirc.load_segment(fn, informat, [0,4096,0,4096], [k], False).astype(numpy.float32)/(ngroups-1)
      # reference pixel subtraction
      for i in range(4096):
        p = numpy.median(numpy.concatenate((persistence_map[f,i,:4],persistence_map[f,i,-4:])))
        persistence_map[f,i,:] = persistence_map[f,i,:] - p
      for i in range(32):
        xmin = 128*i; xmax = xmin+128
        p = numpy.median(persistence_map[f,-4:,xmin:xmax])        
        persistence_map[f,:,xmin:xmax] = persistence_map[f,:,xmin:xmax] - p
    else:
      print('Error: failed to match line {:d}: '.format(f) + content[f])
      exit()
  Q = numpy.asarray(Q)
  # convert DN -> e, subtract dark
  for j in range(NP):
    # gain conversion
    for iy in range(ny):
      ymin = iy*4096//ny; ymax = ymin+4096//ny
      for ix in range(nx):
        xmin = ix*4096//nx; xmax = xmin+4096//nx
        persistence_map[j,ymin:ymax,xmin:xmax] *= gain[iy,ix]
    persistence_map[j,:,:] = persistence_map[j,:,:] - tgroup*dark_current*(ngroups/2.)
    # convert from e in tgroup to e per ln t
    X = numpy.average( numpy.log(1.+tgroup/frtime*numpy.linspace(1,ngroups-1,ngroups-1)) )
    persistence_map[j,:,:] /= X
  print('average ln (tmax/tmin) =', X, ' used ngroups=', ngroups)
  for j in range(NP):
    st = 'pers at {:8.1f} :'.format(Q[j])
    for i in range(1,10): st = st+' {:5.1f}'.format(numpy.percentile(persistence_map[j,:,:], i*10))
    st += ' (unfilt)'
    print(st); this_metadata.append(st)
  #
  # median filter the image
  print('median filtering persistence images ...')
  nfilt = 5
  for j in range(NP): persistence_map[j,4:-4,4:-4] = ndimage.median_filter(persistence_map[j,4:-4,4:-4], size=nfilt, mode='reflect')
  for j in range(NP):
    st = 'pers at {:8.1f} :'.format(Q[j])
    for i in range(1,10): st = st+' {:5.1f}'.format(numpy.percentile(persistence_map[j,:,:], i*10))
    st += ' (filt)'
    print(st); this_metadata.append(st)
  #
  persist_hdu = fits.ImageHDU(persistence_map[1:,:,:])
  persist_hdu.header['EXTNAME'] = 'PERSIST'
  persist_hdu.header['PERSIST'] = (True, 'Persistence implemented')
  for j in range(1,NP):
    persist_hdu.header['Q{:02d}'.format(j)] = (Q[j], 'Stimulus in e')
  alpha = numpy.log(numpy.median(persistence_map[-1,:,:])/numpy.median(persistence_map[-2,:,:]))/numpy.log(Q[-1]/Q[-2])
  print('alpha =', alpha)
  persist_hdu.header['ALPHA'] = (alpha, 'High end exponent')
  persist_hdu.header['SPFILTER'] = (nfilt, 'nxn spatial median')
else:
  this_metadata.append('Skipping persistence ...')
  persist_hdu = fits.ImageHDU(numpy.zeros((2,1,1)))
  persist_hdu.header['EXTNAME'] = 'PERSIST'
  persist_hdu.header['PERSIST'] = (False, 'Does not implement persistence -- placeholder')
  persist_hdu.header['Q01'] = (float(1e4), 'Placeholder')
  persist_hdu.header['Q02'] = (float(1e5), 'Placeholder')
  persist_hdu.header['ALPHA'] = (0., 'Placeholder')
this_metadata.append('')

#######################################################################
# General output
#######################################################################

# Make gain map
gain_hdu= fits.ImageHDU(gain.astype(numpy.float32))
gain_hdu.header['EXTNAME'] = 'GAIN'

print(numpy.shape(badpix), badpix.dtype)
badpix_hdu = fits.ImageHDU(badpix)
badpix_hdu.header['EXTNAME'] = 'BADPIX'
badpix_hdu.header['BIT00'] = 'Disconnected or low response pixel'
badpix_hdu.header['BIT01'] = 'Hot pixel (used TDARK2)'
badpix_hdu.header['BIT02'] = 'Very hot pixel (used TDARK1)'
badpix_hdu.header['BIT03'] = 'Adjacent to pixel with strange response'

# Source information
for k in range(len(this_metadata)): this_metadata[k] = this_metadata[k][:512]
metadata_col = fits.Column(name='INFO', format='512A', array=this_metadata)
src_hdu = fits.BinTableHDU.from_columns([metadata_col])
src_hdu.header['EXTNAME'] = 'SOURCES'

# Final output step
hdul = fits.HDUList([primary_hdu, src_hdu, relqe1_hdu, bfe_hdu, dark_hdu, darkvar_hdu, persist_hdu,\
  saturate_hdu, cnl_hdu, ipc_hdu, ipcflat_hdu, vtpe_hdu, badpix_hdu, noise_hdu, gain_hdu, bias_hdu])
hdul.writeto(configInfo.OUT, overwrite=True)
