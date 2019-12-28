# Routine to process the sequence of flats and darks
# in the WFIRST acceptance testing

import sys
import os, fnmatch
import time
import re
import numpy
import pyirc
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

use_cmap = 'gnuplot'

# Read input information

if len(sys.argv)<4:
  print('Usage: python fdseq.py <input dir> <output stem> <frame numbers> <options>')
  print(len(sys.argv)-1, 'arguments given')
  exit()

input_dir = sys.argv[1]
outstem = sys.argv[2]
m = re.search(r'(\d+),(\d+)', sys.argv[3])
fr1 = int(m.group(1)); fr2 = int(m.group(2))
my_options = sys.argv[4:] # list


# Set parameters

formatpars = 4
nx = ny = 32
subtr_ref = False
sca = 'xxxxx'

# Changes set by options

for this_opt in my_options:

  m = re.search(r'^-n=(\d+),(\d+)$', this_opt)
  if m:
    nx = int(m.group(1))
    ny = int(m.group(2))
    print('Set sizes=',nx,ny)
  m = re.search(r'^-r', this_opt)
  if m: subtr_ref = True
  m = re.search(r'^-s=(\d+)$', this_opt)
  if m: sca = m.group(1)

# Basic calculations

# Size of a block
N = pyirc.get_nside(formatpars)
# Side lengths
dx = N//nx
dy = N//ny

# Now get file list
flist1 = os.listdir(input_dir)
flist = []
slist = []
elist = []
for filename in flist1:
  if fnmatch.fnmatch(filename, '*1400nm*.fits'):
    flist.append(filename)
    s = e = 65535
    m = re.search(r'_(\d+).fits', filename)
    if m: e = int(m.group(1))
    m = re.search(r'ch(\d+)_', filename)
    if m: s = int(m.group(1))
    #
    slist.append(s)
    elist.append(e)
# Sort in order of set & exposure number
#
u = sorted(zip(slist,elist,flist))
Xfiles = [x for _,_,x in u]
Xsets = [x for x,_,_ in u]
Xexposures = [x for _,x,_ in u]
NF = len(Xfiles)
# cleanup stuff we don't need
del u, flist, slist, elist, flist1

print(NF, 'files')
print(Xfiles)
print('Sets ->', Xsets)
print('Exposures ->', Xexposures)
print('')

# Read the reference file
Smed = numpy.zeros((NF,ny,nx))
for k in range(NF):
  print('Reading file',k,'->',Xfiles[k])
  Stot = pyirc.load_segment(input_dir+'/'+Xfiles[k], formatpars, [0,N-1,0,N-1], [fr1,fr2], False)
  S = Stot[0,:,:] - Stot[1,:,:]
  del Stot
  for j in range(ny):
    ref_left = numpy.median(S[j*dy:j*dy+dy,:4])
    ref_right = numpy.median(S[j*dy:j*dy+dy,-4:])
    for i in range(nx):
      Smed[k,j,i] = numpy.median(S[j*dy:j*dy+dy,i*dx:i*dx+dx])
    if subtr_ref:
      Smed[k,j,:] -= (ref_left+ref_right)/2.

#kref=NF-1
#while Xsets[kref]%2==0 or Xexposures[kref]>1: kref-=1
kref = 0
while Xsets[kref]%2==0: kref+=1
print('Ref. exposure = {:3d} (S{:2d}, E{:2d})'.format(kref, Xsets[kref], Xexposures[kref]))

print('Blocks of medians ->')
print(Smed[kref,:,:])
print('')

# Generate ratios of flats and percentiles
plist = [0., 2.28, 15.87, 50.00, 84.13, 97.72, 100.]
np = len(plist)
R = numpy.zeros((NF,ny,nx))
Rptiles = numpy.zeros((NF,np))
for k in range(NF):
  R[k,:,:] = numpy.where(Smed[kref,:,:]>0, Smed[k,:,:]/Smed[kref,:,:], 0*Smed[kref,:,:]-1)
  for ip in range(np):
    temp_array = numpy.where(Smed[kref,:,:]>0, Smed[k,:,:]/Smed[kref,:,:], 0*Smed[kref,:,:]+numpy.nan)
    Rptiles[k,ip] = numpy.nanpercentile(temp_array, plist[ip])

print('Ratios:')
for k in range(NF):
  print('{:3d} {:3d}'.format(Xsets[k], Xexposures[k]), end='')
  for ip in range(np):
    print(' {:8.5f}'.format(Rptiles[k,ip]), end='')
  print('')

# Generate output tables
thisOut = open(outstem+'_ratios.txt', 'w')
thisOut.write('# File: '+outstem+'_ratios.txt\n')
thisOut.write('# This summary created at {:s}\n'.format(time.asctime(time.localtime(time.time()))))
thisOut.write('# Frame numbers {:2d} {:2d}\n'.format(fr1,fr2))
thisOut.write('# Options:')
for x in my_options:
  thisOut.write(' '+x)
thisOut.write('\n')
thisOut.write('# Format: <exposure id> <set> <exposure> <percentiles of R> (7 columns) <R(0,0)> <R(0,1)> ... <R(ny-1,nx-1)>\n')
thisOut.write('#\n')
for k in range(NF):
  thisOut.write('{:3d} {:2d} {:2d}'.format(k, Xsets[k], Xexposures[k]))
  for ip in range(np):
    thisOut.write(' {:8.5f}'.format(Rptiles[k,ip]))
  for iy in range(ny):
    for ix in range(nx):
      thisOut.write(' {:8.5f}'.format(R[k,iy,ix]))
  thisOut.write('\n')
thisOut.close()

# FITS cube output
hdu = fits.PrimaryHDU(R)
hdu.header['DATE'] = format(time.asctime(time.localtime(time.time())))
hdu.header['SCA'] = sca
hdu.header['ORIGIN'] = 'fdseq.py'
hdu.header['FILENAME'] = outstem+'_Rcube.fits'
hdul = fits.HDUList([hdu])
hdul.writeto(outstem+'_Rcube.fits', overwrite=True)

# Generate map of the 2nd frame effect
R2nd = numpy.zeros((ny,nx))
use2nd = 0
for k in range(NF):
  if Xsets[k]%2==1 and Xexposures[k]==2:
    R2nd += R[k,:,:]-R[k-1,:,:]
    use2nd += 1
R2nd /= use2nd
print('Making map of 2nd frame effect with {:d} 2nd frames'.format(use2nd))
print(R2nd)
ar = nx/(ny+0.0)
matplotlib.rcParams.update({'font.size': 8})
F = plt.figure(figsize=(4,3))
S = F.add_subplot(1,1,1)
S.set_title(r'2nd frame signal (%), SCA{:s}'.format(sca))
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(R2nd*100, cmap=use_cmap, aspect=ar, interpolation='nearest', origin='lower',
  vmin=-.7, vmax=.7)
F.colorbar(im, orientation='vertical')
F.set_tight_layout(True)
F.savefig(outstem+'_2nd.pdf')
plt.close(F)

# Generate map of the secular effect
kref1 = 0
while Xsets[kref1]%2==0: kref1+=1
kref2 = NF-1
while Xsets[kref2]%2==0 or Xexposures[kref2]>1: kref2-=1
print('Making map of secular effect')
ar = nx/(ny+0.0)
matplotlib.rcParams.update({'font.size': 8})
F = plt.figure(figsize=(4,3))
S = F.add_subplot(1,1,1)
S.set_title(r'Secular drift, SCA{:s}, Set {:d} -> {:d} (%)'.format(sca,Xsets[kref1],Xsets[kref2]))
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(numpy.log(R[kref2,:,:]/R[kref1,:,:])*100, cmap=use_cmap, aspect=ar, interpolation='nearest', origin='lower',
  vmin=-1, vmax=1)
F.colorbar(im, orientation='vertical')
F.set_tight_layout(True)
F.savefig(outstem+'_sec.pdf')
plt.close(F)

# Generate persistence map
kref1 = 0
while Xsets[kref1]%2==0: kref1+=1
kref2 = kref1
while Xsets[kref2]<2: kref2+=1
print('Making persistence map')
ar = nx/(ny+0.0)
matplotlib.rcParams.update({'font.size': 8})
F = plt.figure(figsize=(4,3))
S = F.add_subplot(1,1,1)
S.set_title(r'Persistence map, SCA{:s}, Set {:d} -> {:d} (%)'.format(sca,Xsets[kref1],Xsets[kref2]))
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(R[kref2,:,:]/R[kref1,:,:]*100, cmap=use_cmap, aspect=ar, interpolation='nearest', origin='lower',
  vmin=-.7, vmax=.7)
F.colorbar(im, orientation='vertical')
F.set_tight_layout(True)
F.savefig(outstem+'_per.pdf')
plt.close(F)

# Response in flat/dark sequences
F = plt.figure(figsize=(6,4))
S = F.add_subplot(1,1,1)
S.set_title('flat/dark signals: SCA' + sca)
S.set_xlim(0,NF)
S.set_ylim(-.02,.02)
C=[]
for k in range(NF):
  if Xexposures[k]==1: C.append(k)
S.tick_params(axis='x', labelrotation=90)
S.xaxis.set_ticks(C); del C
S.set_xlabel('Exposure number')
S.set_ylabel('Relative response')
myX = numpy.array(range(NF))+.5
S.yaxis.set_ticks(numpy.linspace(-.02,.02,num=9))
S.grid(True, color='g', linestyle='-', linewidth=.333)
S.scatter(myX, Rptiles[:,np//2]-1, s=1, marker='D', color='r')
S.scatter(myX, Rptiles[:,0]-1, s=1, marker=10, color='r')
S.scatter(myX, Rptiles[:,-1]-1, s=1, marker=11, color='r')
for ip in range(1,np-1):
  if ip!=np//2:
    for o in range(2):
      S.scatter(myX, Rptiles[:,ip]-o, s=1, marker='.', color='b')
S.scatter(myX, Rptiles[:,np//2], s=1, marker='x', color='k')
S.scatter(myX, Rptiles[:,0], s=1, marker='_', color='k')
S.scatter(myX, Rptiles[:,-1], s=1, marker='_', color='k')
F.set_tight_layout(True)
F.savefig(outstem+'_fdevol.pdf')
plt.close(F)

