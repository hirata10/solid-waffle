import sys
import numpy
import re
import mtfutils
from astropy.io import fits

inputdata = []
nrepeat = []
dir = []
cppinit = []

outstem = ''

nside = 4096

# Get configuration
config_file = sys.argv[1]
with open(config_file) as myf: content = myf.read().splitlines()
is_input = False
for line in content:
  # Cancellations
  m = re.search(r'^[A-Z]+\:', line)
  if m: is_input = False

  # Searches for files
  if is_input:
    m = re.search(r'^\s*(\S+)\s+(\d+)\s+(\w)\s+(\S+)$', line)
    if m:
      inputdata += [m.group(1)]
      nrepeat += [int(m.group(2))]
      dir += [m.group(3).upper()]
      cppinit += [float(m.group(4))]

  m = re.search(r'^INPUT\:', line)
  if m: is_input = True

  # Bin sizes
  m = re.search(r'^NBIN:\s*(\d+)\s+(\d+)', line)
  if m:
    nx = int(m.group(1))
    ny = int(m.group(2))

  m = re.search(r'^OUTPUT:\s*(\S+)', line)
  if m: outstem = m.group(1)

# end read config

if outstem=='':
  print('Error: no output specified.')
  exit()

n_input_group = len(inputdata)
max_nrepeat = numpy.amax(numpy.asarray(nrepeat))
dx = nside//nx; dy = nside//ny

# Tell the user what we learned
print('Using files:')
for k in range(n_input_group): print('   ',inputdata[k], 'repeats:', nrepeat[k])
print('max. number of repeats =', max_nrepeat)
print('number of superpixels: {:d} (in y) x {:d} (in x) = {:d} (total)'.format(ny,nx,nx*ny))

# Loop over the superpixels
all_cpp = numpy.zeros((n_input_group, ny, nx))
all_mtf2 = numpy.zeros((n_input_group, ny, nx))
for iy in range(ny):
  for ix in range(nx):
    # get super pixel boundaries
    xmin = dx*ix; xmax = dx*(ix+1); ymin = dy*iy; ymax = dy*(iy+1)
    if xmin<4: xmin=4
    if xmax>nside-4: xmax=nside-4
    if ymin<4: ymin=4
    if ymax>nside-4: ymax=nside-4
    print('processing super pixel {:d},{:d}'.format(iy,ix), '--> [{:4d}:{:4d},{:4d}:{:4d}]'.format(ymin,ymax,xmin,xmax))
    wx = xmax-xmin; wy = ymax-ymin # widths

    # make array of subimages
    data = numpy.zeros((n_input_group, max_nrepeat, wy, wx))
    for k in range(n_input_group):
      for ik in range(nrepeat[k]):
        hdul = fits.open(re.sub(r'\*', '{:03d}'.format(ik+1), inputdata[k])) 
        if len(hdul)==1:
          data[k,ik,:,:] = hdul[0].data[ymin:ymax,xmin:xmax]
        if len(hdul)==2:
          data[k,ik,:,:] = hdul[1].data[0,0,ymin:ymax,xmin:xmax].astype(numpy.float64) - hdul[1].data[0,-1,ymin:ymax,xmin:xmax].astype(numpy.float64)
        hdul.close()

    # power spectrum analysis
    PSH = numpy.zeros((n_input_group,wx))
    PSV = numpy.zeros((n_input_group,wy))
    for k in range(n_input_group):
      for ik in range(nrepeat[k]):
        PS = numpy.abs(numpy.fft.ifft2(data[k,ik,:,:]))**2
        PS[0,:] = 0.; PS[:,0] = 0.
        PSH[k,:] = PSH[k,:] + numpy.sum(PS, axis=0)
        PSV[k,:] = PSV[k,:] + numpy.sum(PS, axis=1)
      PSuse = numpy.copy(PSH[k,:])
      w = wx
      if dir[k]=='V':
        PSuse = numpy.copy(PSV[k,:])
        w = wy
      if iy == 0 and ix == 0 and k == 0 :
        all_PSuse = numpy.zeros((ny,nx,n_input_group,numpy.size(PSuse)))
      all_PSuse[iy,ix,k,:] = PSuse
        
      c0,cw,a1,a2,res = mtfutils.get_triangle_from_ps(PSuse,cppinit[k])
      print('{:2d} {:3d} {:3d} {:6.4f} {:6.4f} {:10.4E} {:7.5f} {:9.4E}'.format(k, iy, ix, c0,cw,a1,a2/a1,res))
      # put in cube
      all_cpp[k,iy,ix] = c0
      all_mtf2[k,iy,ix] = a2/a1
      
all_PSuse.reshape((ny*nx, n_input_group, numpy.size(PSuse)))
hdu = fits.PrimaryHDU(all_PSuse)
hdulist = fits.HDUList([hdu])
hdulist.writeto('PS_speckle.fits', overwrite = True)

# make output picture
# 2 slices -- wavenumber & MTF**2
nslice = 2
out_array = numpy.zeros((nslice,n_input_group,ny,nx))
for k in range(n_input_group):
  out_array[0,k,:,:] = all_cpp[k,:,:]
  out_array[1,k,:,:] = all_mtf2[k,:,:]
hdu = fits.PrimaryHDU(numpy.transpose(out_array, axes=(0,2,1,3)).reshape((nslice*ny,n_input_group*nx)))
hdu.writeto(outstem+'info.fits', overwrite=True)

# report summary
print('== Summary of results ==')
print('k [cpp]|  mtf  |std mtf')
for k in range(n_input_group):
  m2_50 = numpy.median(all_mtf2[k,:,:])
  m2_25 = numpy.percentile(all_mtf2[k,:,:],25)
  m2_75 = numpy.percentile(all_mtf2[k,:,:],75)
  print('{:7.5f} {:7.5f} {:7.5f}'.format(numpy.mean(all_cpp[k,:,:]),
    numpy.sqrt(m2_50), (numpy.sqrt(m2_75)-numpy.sqrt(m2_25))/1.349/numpy.sqrt(numpy.size(all_mtf2[k,:,:])-1.)*numpy.sqrt(numpy.pi/2.)))
