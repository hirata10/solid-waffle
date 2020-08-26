import numpy
from astropy.io import fits
import sys
import re
import pyirc
import time

thisversion = 4

Narg = len(sys.argv)
if Narg<3:
  print('Error: {:d} arguments\n'.format(Narg))
  print('Usage: python spr_reduce.py <input file> <output stem> <options>\n')
  print('Options:')
  print('  -f=<format> -> controls input format')
  print('  -n=<nfile> -> controls number of files used (assumed sequential if n>=2)')
  print('  -p=<pattern> -> controls pattern number')
  print('  -d=<dark file> -> dark data cube, same format, used for masking')
  print('  -nd=<ndarks> -> controls number of dark files used (assumed sequential if n>=2)')
  print('  -sca=<scanum> -> SCA number (for output)')
  print('  -sd -> subtract dark (not recommended at this time, hasn\'t worked as well as we had hoped)')
  print('  -i -> interpolate masked pixels')
  print('  -a=<parameter> -> alternate file naming scheme')
  print('  -nl=<summary file> -> summary file for NL information')
  exit()

nfile = 1
formatpars = 4
tdark = 10
outstem = sys.argv[2]
ipc_pattern = 1
usedark = False; ndark = 1; darkfile = ''
subtr_dark = False
interp_alpha = False
sca = 'xxxxx'
name_scheme = 1
use_nl = False

for k in range(3,Narg):
  m = re.search(r'^-f=(\d+)$', sys.argv[k])
  if m: formatpars = int(m.group(1))
  m = re.search(r'^-n=(\d+)$', sys.argv[k])
  if m: nfile = int(m.group(1))
  m = re.search(r'^-p=(\d+)$', sys.argv[k])
  if m: ipc_pattern = int(m.group(1))
  m = re.search(r'^-d=(\S+)$', sys.argv[k])
  if m:
    usedark = True
    darkfile = m.group(1)
  m = re.search(r'^-nd=(\d+)$', sys.argv[k])
  if m: ndark = int(m.group(1))
  m = re.search(r'^-sd$', sys.argv[k])
  if m: subtr_dark = True
  m = re.search(r'^-i$', sys.argv[k])
  if m: interp_alpha = True
  m = re.search(r'^-sca=(\d+)', sys.argv[k])
  if m: sca = m.group(1)
  m = re.search(r'^-a=(\d+)', sys.argv[k])
  if m: name_scheme = int(m.group(1))
  m = re.search(r'^-nl=(\S+)$', sys.argv[k])
  if m:
    use_nl = True
    nlfile = m.group(1)

N = pyirc.get_nside(formatpars)
dmap = numpy.zeros((nfile,N,N))

# Pull down information from NL file
if use_nl:
  summaryinfo = numpy.loadtxt(nlfile)
  sum_nx = int(numpy.amax(summaryinfo[:,0]))+1
  sum_ny = int(numpy.amax(summaryinfo[:,1]))+1
  omax = 6; colindex = [0]*(omax+1)
  nlcoefs = numpy.zeros((omax+1, sum_ny, sum_nx))
  f = open(nlfile, 'r')
  for x in f:
    m = re.search('^\# +(\d+), additional non-linearity coefficient, order (\d+) ', x)
    if m:
      thiscol = int(m.group(1))
      thisord = int(m.group(2))
      colindex[thisord] = thiscol
      if thisord>omax:
        print('Error: increase omax={:d} to {:d}'.format(omax,thisord))
        exit()
  print('NL information ->', sum_ny, sum_nx, omax, colindex)
  for p in range(1,omax+1):
    if colindex[p]>0:
      nlcoefs[p,:,:] = summaryinfo[:,colindex[p]-1].reshape(sum_ny,sum_nx)
  f.close()

# IPC patterns:
#
dx = dy = 8
if ipc_pattern==1:
  dx = dy = 8
# < alternate dx,dy would go here >

nx = N//dx; ny = N//dy
rx = numpy.zeros(nx, dtype=int)
ry = numpy.zeros(ny, dtype=int)

signals = numpy.zeros((nfile,13,ny,nx))
dmask = numpy.zeros((ny,nx))

# List of IPC patterns (actually a list)
if ipc_pattern==1:
  for j in range(16):
    rx[32*j:32*j+16] = numpy.arange(0,128,8) + 256*j + 8
    rx[32*j+16:32*j+32] = 256*j + 247 - numpy.arange(0,128,8)[::-1]
  ry = numpy.arange(dy-1,N,dy)

# Make dark map
if usedark:
  # Dark map
  D = numpy.zeros((ndark,N,N))
  for j in range(ndark):
    thisfile = darkfile + ''
    if j>0:
      m = re.search(r'(.+_)(\d+)\.fits$', darkfile)
      if m:
        new_index = int(m.group(2)) + j
        thisfile = m.group(1) + '{:03d}.fits'.format(new_index)
      else:
        print('Error: can\'t construct new file name.')
        exit()
    thisframe = pyirc.load_segment(thisfile, formatpars, [0,N,0,N], [1,1+tdark], True)
    D[j,:,:] = thisframe[0,:,:] - thisframe[1,:,:]
  darkframe = numpy.median(D,axis=0)/float(tdark)
  del D

  # Make FITS output of dark map
  hdu = fits.PrimaryHDU(darkframe)
  hdu.header['DATE'] = format(time.asctime(time.localtime(time.time())))
  hdu.header['SCA'] = sca
  hdu.header['ORIGIN'] = 'spr_reduce.py'
  hdu.header['VERSION'] = '{:d}'.format(thisversion)
  hdu.header['FILETYPE'] = 'dark map for IPC masking'
  hdul = fits.HDUList([hdu])
  hdul.writeto(outstem + '_sprdark.fits', overwrite=True)

filelist = []
for j in range(nfile):
  thisfile = sys.argv[1]
  if j>0:
    if name_scheme==1:
      m = re.search(r'(.+_)(\d+)\.fits$', sys.argv[1])
      if m:
        new_index = int(m.group(2)) + j
        thisfile = m.group(1) + '{:03d}.fits'.format(new_index)
      else:
        print('Error: can\'t construct new file name.')
        exit()
    if name_scheme==2:
      m = re.search(r'(.+)_1_(.+\.fits)$', sys.argv[1])
      if m:
        thisfile = m.group(1) + '_{:d}_'.format(j+1) + m.group(2)
      else:
        print('Error: can\'t construct new file name.')
        exit()

  filelist.append(thisfile)
  thisframe = pyirc.load_segment(thisfile, formatpars, [0,N,0,N], [1,2], True)
  dmap[j,:,:] = thisframe[0,:,:] - thisframe[1,:,:]
  if subtr_dark: dmap[j,:,:] -= darkframe

  # Non-linearity correction, if used
  if use_nl:
    sys.stdout.write('Non-linearity corrections: super-rows (of {:d}): '.format(sum_ny)); sys.stdout.flush()
    for iys in range(sum_ny):
      sys.stdout.write('{:d} '.format(iys)); sys.stdout.flush()
      ysmin = iys*(N//sum_ny)
      ysmax = (iys+1)*(N//sum_ny)
      for ixs in range(sum_nx):
        xsmin = ixs*(N//sum_nx)
        xsmax = (ixs+1)*(N//sum_nx)
        S = dmap[j,ysmin:ysmax,xsmin:xsmax] # makes subarray
        Sf = numpy.copy(S)
        # want to solve Sf = S + c_2 S^2 + c_3 S^3 + ...
        for k in range(20):
          # iterative solution
          Sp = numpy.copy(S)
          for p in range(2,omax+1):
            Sp += nlcoefs[p,iys,ixs]*S**p
          S += Sf-Sp
    print('Done.')

  # subtractions to get "signals" (background, center, horiz, vert, diag)
  for iy in range(ny):
    yc = ry[iy]
    if yc>=5 and yc<N-5:
      for ix in range(nx):
        xc = rx[ix]
        if xc>=5 and xc<N-5:
          signals[j,0,iy,ix] = numpy.median(numpy.concatenate(( dmap[j,yc-3:yc-1,xc-3:xc+4].flatten(), dmap[j,yc+2:yc+4,xc-3:xc+4].flatten(),
                                 dmap[j,yc-1:yc+2,xc-3:xc-1].flatten(), dmap[j,yc-1:yc+2,xc+2:xc+4].flatten() )))
          signals[j,1,iy,ix] = dmap[j,yc,xc] - signals[j,0,iy,ix]
          # horizontal
          signals[j,2,iy,ix] = (dmap[j,yc,xc-1]+dmap[j,yc,xc+1])/2. - signals[j,0,iy,ix]
          # vertical
          signals[j,3,iy,ix] = (dmap[j,yc-1,xc]+dmap[j,yc+1,xc])/2. - signals[j,0,iy,ix]
          # diagonal
          signals[j,4,iy,ix] = (dmap[j,yc-1,xc-1]+dmap[j,yc+1,xc-1]+dmap[j,yc-1,xc+1]+dmap[j,yc+1,xc+1])/4. - signals[j,0,iy,ix]
          # individual pixels, going from "right" to "upper-right" and then around
          signals[j,5,iy,ix] = dmap[j,yc,xc+1] - signals[j,0,iy,ix]
          signals[j,6,iy,ix] = dmap[j,yc+1,xc+1] - signals[j,0,iy,ix]
          signals[j,7,iy,ix] = dmap[j,yc+1,xc] - signals[j,0,iy,ix]
          signals[j,8,iy,ix] = dmap[j,yc+1,xc-1] - signals[j,0,iy,ix]
          signals[j,9,iy,ix] = dmap[j,yc,xc-1] - signals[j,0,iy,ix]
          signals[j,10,iy,ix] = dmap[j,yc-1,xc-1] - signals[j,0,iy,ix]
          signals[j,11,iy,ix] = dmap[j,yc-1,xc] - signals[j,0,iy,ix]
          signals[j,12,iy,ix] = dmap[j,yc-1,xc+1] - signals[j,0,iy,ix]

# Make FITS output of difference map
hdu = fits.PrimaryHDU(numpy.mean(dmap,axis=0))
hdu.header['DATE'] = format(time.asctime(time.localtime(time.time())))
hdu.header['SCA'] = sca
hdu.header['ORIGIN'] = 'spr_reduce.py'
hdu.header['VERSION'] = '{:d}'.format(thisversion)
hdu.header['FILETYPE'] = 'SPR difference map'
hdul = fits.HDUList([hdu])
hdul.writeto(outstem + '_sprmean.fits', overwrite=True)

medsignals = numpy.median(signals,axis=0)
print('median signal = {:8.6f} DN'.format(numpy.median(medsignals[1,:,:])))
alpha = numpy.zeros((13,ny,nx))

# Masking based on the dark
if usedark:
  for iy in range(ny):
    yc = ry[iy]
    if yc>=5 and yc<N-5:
      for ix in range(nx):
        xc = rx[ix]
        if xc>=5 and xc<N-5:
          if (numpy.amax(darkframe[yc-1:yc+2,xc-1:xc+2])>1e-3*tdark*medsignals[1,iy,ix]):
            dmask[iy,ix] = 1

# mask if central pixel isn't the brightest by at least a factor of 10
# -- would usually indicate a problem
for iy in range(ny):
  for ix in range(nx):
    if .1*medsignals[1,iy,ix]<numpy.amax(medsignals[5:13,iy,ix]):
      dmask[iy,ix] = 1

# alpha map
den = medsignals[1,:,:] + 2*medsignals[2,:,:] + 2*medsignals[3,:,:] + 4*medsignals[4,:,:] + 1e-99
# (the 1e-99 ensures that 0's are passed through)
alpha[:11,:,:] = medsignals[2:,:,:]/den
alpha[11,:,:] = (alpha[0,:,:] + alpha[1,:,:])/2.
alpha[12,:,:] = dmask[:,:]

# filling in
if ipc_pattern==1:
  alpha[:,-1,:] = alpha[:,-2,:]
  for j in range(16):
    alpha[:,:,15+32*j] = alpha[:,:,14+32*j]
    alpha[:,:,16+32*j] = alpha[:,:,17+32*j]

# interpolation
for iy in range(ny):
  for ix in range(nx):
    if dmask[iy,ix]>.5:
      # first try 8 nearest neighbors
      aDen = numpy.sum(1-dmask[iy-1:iy+2,ix-1:ix+2])
      if aDen>0:
        for k in range(12):
          alpha[k,iy,ix] = numpy.sum((1-dmask[iy-1:iy+2,ix-1:ix+2])*alpha[k,iy-1:iy+2,ix-1:ix+2])/aDen

# Make FITS output of IPC map
hdu = fits.PrimaryHDU(alpha)
hdu.header['DATE'] = format(time.asctime(time.localtime(time.time())))
hdu.header['SCA'] = sca
hdu.header['ORIGIN'] = 'spr_reduce.py'
hdu.header['VERSION'] = '{:d}'.format(thisversion)
hdu.header['FILETYPE'] = 'IPC data cube'
hdu.header['DX'] = str(dx)
hdu.header['DY'] = str(dy)
hdu.header['SLICE01'] = 'alpha_horizontal'
hdu.header['SLICE02'] = 'alpha_vertical'
hdu.header['SLICE03'] = 'alpha_diagonal'
hdu.header['SLICE04'] = 'alpha dx=+1, dy= 0'
hdu.header['SLICE05'] = 'alpha dx=+1, dy=+1'
hdu.header['SLICE06'] = 'alpha dx= 0, dy=+1'
hdu.header['SLICE07'] = 'alpha dx=-1, dy=+1'
hdu.header['SLICE08'] = 'alpha dx=-1, dy= 0'
hdu.header['SLICE09'] = 'alpha dx=-1, dy=-1'
hdu.header['SLICE10'] = 'alpha dx= 0, dy=-1'
hdu.header['SLICE11'] = 'alpha dx=+1, dy=-1'
hdu.header['SLICE12'] = 'alpha (average of 4 nearest neighbors)'
hdu.header['SLICE13'] = 'mask (0 = normal, 1 = masked), not foolproof'
for k in range(Narg):
  keyword = 'ARGV{:02d}'.format(k)
  hdu.header[keyword] = sys.argv[k]
hdu.header['MASKSIZE'] = '{:d}/{:d}'.format(int(numpy.sum(dmask)), nx*ny)
hdu.header['MEDSIG'] = ('{:8.2f}'.format(numpy.median(medsignals[1,:,:])), 'Median signal in central pixel')
for k in range(len(filelist)):
  keyword = 'INF{:02d}'.format(k)
  hdu.header[keyword] = filelist[k]
hdul = fits.HDUList([hdu])
hdul.writeto(outstem + '_alpha.fits', overwrite=True)

print('median alpha information ->')
print(numpy.median(alpha,axis=[1,2]))

print('Number of masked pixels = {:d}/{:d}'.format(int(numpy.sum(dmask)), nx*ny))
