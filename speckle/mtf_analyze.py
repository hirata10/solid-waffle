import sys
import numpy
import numpy.linalg
import re
import mtfutils
from astropy.io import fits

inputdata = []
nrepeat = []
dir = []
cppinit = []
wavelen = []
d = []
s = []
z = 127
P = 0.01

outstem = ''

nskip = 1 # skip beginning & end of rows in power spectrum before stacking

finalframe = -1 # last frame to use

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
    m = re.search(r'^\s*(\S+)\s+(\d+)\s+(\w)\s+(\S+)\s+(\S+)\s+(\d+\.\d+)\s+(\S+)$', line)
    if m:
      inputdata += [m.group(1)]
      nrepeat += [int(m.group(2))]
      dir += [m.group(3).upper()]
      cppinit += [float(m.group(4))]
      wavelen += [float(m.group(5))] # wavelength
      d += [float(m.group(6))] # slit spacing
      s += [float(m.group(7))] # slit width
  
  m = re.search(r'^INPUT\:', line)
  if m: is_input = True

  # Bin sizes
  m = re.search(r'^NBIN:\s*(\d+)\s+(\d+)', line)
  if m:
    nx = int(m.group(1))
    ny = int(m.group(2))

  # Bin sizes
  m = re.search(r'^NSKIP:\s*(\d+)', line)
  if m:
    nskip = int(m.group(1))

  # Bin sizes
  m = re.search(r'^FINALFRAME:\s*(\d+)', line) # input in DS9 ordering
  if m:
    finalframe = int(m.group(1))-1 # convert to Python

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

# create empty Mxx matrix
Mxx = numpy.zeros((ny,nx))

# loop over the superpixels
all_cpp = numpy.zeros((n_input_group, ny, nx))
all_mtf2 = numpy.zeros((n_input_group, ny, nx))
all_cpp_force = numpy.zeros((n_input_group, ny, nx))
all_dcpp_force = numpy.zeros((n_input_group, ny, nx))
all_mtf2_force = numpy.zeros((n_input_group, ny, nx))
all_mtf2_force_subregion = numpy.zeros((n_input_group, ny, nx))
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
          data[k,ik,:,:] = hdul[1].data[0,0,ymin:ymax,xmin:xmax].astype(numpy.float64) - hdul[1].data[0,finalframe,ymin:ymax,xmin:xmax].astype(numpy.float64)
        hdul.close()

    # create values for Mxx
    x = (wx/2) + xmin
    x = (x-(nside/2))*P
    y = (wy/2) + ymin
    y = (y-(nside/2))*P
    M = ((1 + ((x**2 + y**2)/z**2))**(-3/2))*(1 + ((y**2/z**2)))
    Mxx[iy,ix] = M

    
    # power spectrum analysis
    PSH = numpy.zeros((n_input_group,wx))
    PSV = numpy.zeros((n_input_group,wy))
    for k in range(n_input_group):
      # calculate ideal u and delta u
      u_ideal = ((d[k]*P)/(z*wavelen[k]))*M
      delta_u_ideal = ((s[k]*P)/(z*wavelen[k]))*M
      for ik in range(nrepeat[k]):
        PS = numpy.abs(numpy.fft.ifft2(data[k,ik,:,:]))**2
        PS[0,:] = 0.; PS[:,0] = 0.
        if nskip>1:
          PSH[k,:] = PSH[k,:] + numpy.sum(PS[nskip:-nskip,:], axis=0)
          PSV[k,:] = PSV[k,:] + numpy.sum(PS[:,nskip:-nskip], axis=1)
        else:
          PSH[k,:] = PSH[k,:] + numpy.sum(PS, axis=0)
          PSV[k,:] = PSV[k,:] + numpy.sum(PS, axis=1)
      PSphi = numpy.copy(PSH[k,:]) # multiple copies for various purposes
      PSdis = numpy.copy(PSH[k,:])
      PSuse = numpy.copy(PSH[k,:])
      w = wx
      # get the length of power spectrum without reference pixels
      expected_PSlen = nside//nx
      if dir[k]=='V':
        PSuse = numpy.copy(PSV[k,:])
        PSdis = numpy.copy(PSV[k,:])
        PSphi = numpy.copy(PSV[k,:])
        w = wy
        expected_PSlen = nside//ny
      if iy == 0 and ix == 0 and k == 0 : # executes first so all_PSuse is defined with correct length, and always has reference pixels despite direction of fringes
        all_PSuse = numpy.zeros((ny,nx,n_input_group,numpy.size(PSuse) + 4)) # definition starts in a super pixel with reference pixels, add 4 to use proper length
        PS_display = numpy.zeros_like(all_PSuse)
      if numpy.size(PSuse) < expected_PSlen - 1: # if it does have reference pixels
        PSdis = numpy.interp(numpy.linspace(0,numpy.size(PSuse)-1,expected_PSlen),numpy.linspace(0,numpy.size(PSuse)-1,numpy.size(PSuse)),PSuse)
        #if ix == 0 or ix == 7: # accomodates for both right and left reference pixels
        PSuse = numpy.append(PSuse, (0,0,0,0))

      # add each super pixel power spectrum to overall matrix 
      all_PSuse[iy,ix,k,:] = PSuse
      PS_display[iy,ix,k,:] = PSdis

      # create and (possibly) save the phi matrix for each super pixel
      phi = mtfutils.make_phi_matrix(PSphi,u_ideal,delta_u_ideal)
      phi2 = mtfutils.make_phi_matrix(PSphi,u_ideal,delta_u_ideal,subregion=True)
      hdu = fits.PrimaryHDU(phi)
      hdulist = fits.HDUList([hdu])
      #hdulist.writeto(outstem+'_r{:02d}_y{:02d}_x{:02d}_phi_matrix.fits'.format(k,iy,ix), overwrite = True) 
      hdu = fits.PrimaryHDU(phi2)
      hdulist = fits.HDUList([hdu])
      #hdulist.writeto(outstem+'PHI2_r{:02d}_y{:02d}_x{:02d}_phi_matrix.fits'.format(k,iy,ix), overwrite = True) 
      # ... and get the fit amplitudes
      amplitudes = numpy.linalg.solve(phi@phi.T,phi@PSphi)
      amplitudes2 = numpy.linalg.solve(phi2@phi2.T,phi2@PSphi)
      all_mtf2_force[k,iy,ix] = amplitudes[1]/amplitudes[0]*2
      all_mtf2_force_subregion[k,iy,ix] = amplitudes2[0]/amplitudes[0]*2
      all_cpp_force[k,iy,ix] = u_ideal
      all_dcpp_force[k,iy,ix] = delta_u_ideal
      
      c0,cw,a1,a2,res = mtfutils.get_triangle_from_ps(PSuse,cppinit[k])
      print('{:2d} {:3d} {:3d} {:6.4f} {:6.4f} {:10.4E} {:7.5f} {:9.4E}    {:6.4f} {:10.4E} {:10.4E} {:10.4E} {:10.4E}   {:10.4E} {:10.4E} {:10.4E} {:10.4E}'.format(k, iy, ix, c0,cw,a1,a2/a1,res,
        u_ideal, amplitudes[0], amplitudes[1], amplitudes[2], amplitudes[3], amplitudes2[0], amplitudes2[1], amplitudes2[2], amplitudes2[3]))
      # put in cube
      all_cpp[k,iy,ix] = c0
      all_mtf2[k,iy,ix] = a2/a1

# reshape power spectrum
all_PSuse = all_PSuse.reshape((ny*nx, n_input_group, numpy.size(PSuse)))
all_PSuse = all_PSuse.transpose((1,0,2))
PS_display = PS_display.reshape((ny*nx, n_input_group, numpy.size(PSdis)))
PS_display = PS_display.transpose((1,0,2))

# save power spectrum for analysis
hdu = fits.PrimaryHDU(all_PSuse) 
hdulist = fits.HDUList([hdu]) 
hdulist.writeto(outstem+'_PS_speckle.fits', overwrite = True)

# save power spectrum for display
hdu = fits.PrimaryHDU(PS_display) 
hdulist = fits.HDUList([hdu]) 
hdulist.writeto(outstem+'_PS_speckle_display.fits', overwrite = True)

# save Mxx matrix
hdu = fits.PrimaryHDU(Mxx)
hdulist = fits.HDUList([hdu])
hdulist.writeto('Mxx_ref.fits', overwrite = True) # does not need outstem as it is dependent on the super pixels themselves and not the data

# make output picture
# 2 slices -- wavenumber & MTF**2
nslice = 2
out_array = numpy.zeros((nslice,n_input_group,ny,nx))
for k in range(n_input_group):
  out_array[0,k,:,:] = all_cpp[k,:,:]
  out_array[1,k,:,:] = all_mtf2[k,:,:]
hdu = fits.PrimaryHDU(numpy.transpose(out_array, axes=(0,2,1,3)).reshape((nslice*ny,n_input_group*nx)))
out_array2 = numpy.zeros((nslice,n_input_group,ny,nx))
for k in range(n_input_group):
  out_array2[0,k,:,:] = all_cpp_force[k,:,:]
  out_array2[1,k,:,:] = all_mtf2_force_subregion[k,:,:]
hdu2 = fits.ImageHDU(numpy.transpose(out_array2, axes=(0,2,1,3)).reshape((nslice*ny,n_input_group*nx)))
hdulist = fits.HDUList([hdu,hdu2])
hdulist.writeto(outstem+'info.fits', overwrite=True)

# report summary
print('== Summary of results (all fit) ==')
print('k [cpp]|  mtf  |std mtf')
for k in range(n_input_group):
  m2_50 = numpy.median(all_mtf2[k,:,:])
  m2_25 = numpy.percentile(all_mtf2[k,:,:],25)
  m2_75 = numpy.percentile(all_mtf2[k,:,:],75)
  print('{:7.5f} {:7.5f} {:7.5f}'.format(numpy.mean(all_cpp[k,:,:]),
    numpy.sqrt(m2_50), (numpy.sqrt(m2_75)-numpy.sqrt(m2_25))/1.349/numpy.sqrt(numpy.size(all_mtf2[k,:,:])-1.)*numpy.sqrt(numpy.pi/2.)))

print('')

# report summary of force photometry results
print('== Summary of results (forced fit) ==')
print('k [cpp]|dk[cpp]|  mtf2  |stdmtf2')
for k in range(n_input_group):
  m2_50 = numpy.median(all_mtf2_force[k,:,:])
  m2_25 = numpy.percentile(all_mtf2_force[k,:,:],25)
  m2_75 = numpy.percentile(all_mtf2_force[k,:,:],75)
  print('{:7.5f} {:7.5f} {:8.5f} {:7.5f}'.format(numpy.mean(all_cpp_force[k,:,:]), numpy.mean(all_dcpp_force[k,:,:]),
    m2_50, (m2_75-m2_25)/1.349/numpy.sqrt(numpy.size(all_mtf2[k,:,:])-1.)*numpy.sqrt(numpy.pi/2.)))

print('')

# report summary of force photometry results
print('== Summary of results (forced fit - subregion) ==')
print('k [cpp]|dk[cpp]|  mtf2  |stdmtf2')
for k in range(n_input_group):
  m2_50 = numpy.median(all_mtf2_force_subregion[k,:,:])
  m2_25 = numpy.percentile(all_mtf2_force_subregion[k,:,:],25)
  m2_75 = numpy.percentile(all_mtf2_force_subregion[k,:,:],75)
  print('{:7.5f} {:7.5f} {:8.5f} {:7.5f}'.format(numpy.mean(all_cpp_force[k,:,:]), numpy.mean(all_dcpp_force[k,:,:]),
    m2_50, (m2_75-m2_25)/1.349/numpy.sqrt(numpy.size(all_mtf2[k,:,:])-1.)*numpy.sqrt(numpy.pi/2.)))

print('')
