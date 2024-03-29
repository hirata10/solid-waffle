""" Start with a realization of a perfect detector plus some BFE effect to
the charge in a given pixel
1.  The mean charge <Q_a(i,j)> = It_a + 0.5*Sigma_a*I^2*t_a^2 (Eqn 36) 
2.  Realization of a 4096 x 4096 pixel^2 grid with 66 time samples
This final result data cube will then be [4k, 4k, 66] in dimensions
and will require some kind of identifying header, currently taken
from one of the DCL flats.  This is charge.
3.  Convert to ADU/DN via gain ~ number(e-)/counts;  set ~ FWD/2^16

Run tests with 128x128 or something small...
current I is per pixel (units: e/s)
time t (units: s)
time step for read from t_a to t_a+1 (will need to check # for convergence)

NOTE: to run this, one needs: 
 ** a copy of a DCL flat file with the name set below
"""
import sys
import numpy as np
from numpy.random import randn,poisson,normal
import scipy.signal as signal
import astropy.io.fits as fits
import fitsio
from fitsio import FITS,FITSHDR
import re
sys.path.insert(0, '../')
#sys.path.insert(0, '/users/PCON0003/cond0080/src/solid-waffle/')
from pyirc import *
from ftsolve import p2kernel
from detector_functions import *

# Defaults
formatpars = 1
tsamp = 66
substep = 2
I = 10.0
QE = 0.8
delta_tsamp = 3.0 # arbitrary for now (s)
gain = 1.5 # arbitrary scalar e-/DN
outfile = 'DefaultOutput.fits'
wavemode = 'ir' # options are 'ir' or 'vis'; used only for BFE kernel choice
rngseed = 1000
noisemode = 'none'
bfemode = 'true'
lipcmode = 'false'
lipc_alpha = [0.01]
nlmode = 'false'
nlbeta = 1.5 # (ppm/e-)
# nlcoeffs_arr are [c_2,c_3,c_4] with c_j in units of ppm electrons^(1-j) 
nlcoeffs_arr = [-1.5725,1.9307e-5,-1.4099e-10]
reset_frames = [0]

# quantum yield defaults
QY_omega = 0.
QY_offset = 2
QY_p2 = numpy.zeros((5,5)); QY_p2[2,2] = 1

# Read in information
config_file = sys.argv[1]
with open(config_file) as myf: content = myf.read().splitlines()
for line in content:
  # Format
  m = re.search(r'^FORMAT:\s*(\d+)', line)
  if m: formatpars = int(m.group(1))

  # Number of reads
  m = re.search(r'^NREADS:\s*(\d+)', line)
  if m: tsamp = int(m.group(1))
  # substeps
  m = re.search(r'^SUBSTEPS:\s*(\d+)', line)
  if m: substep = int(m.group(1))

  # Time step (s)
  m = re.search(r'^DT:\s*(\S+)', line)
  if m: delta_tsamp = float(m.group(1))

  # Gain (e/DN)
  m = re.search(r'^GAIN:\s*(\S+)', line)
  if m: gain = float(m.group(1))

  # Quantum yield
  m = re.search(r'^QY:\s*(\d.*)', line)
  if m:
    QY_pars_str = m.group(1).split(" ")
    QY_omega = float(QY_pars_str[0])
    #sig = float(m.group(2))
    QY_cov = [ float(QY_pars_str[x]) for x in range(1, len(QY_pars_str)) ]
    #QY_cxx = float(m.group(2))
    #QY_cxy = float(m.group(3))
    #QY_cyy = float(m.group(4))
    QY_offset = 2
    #QY_p2 = p2kernel(np.asarray([sig**2,0,sig**2]), QY_offset, 256)
    QY_p2 = p2kernel(np.asarray(QY_cov), QY_offset, 256)
  # Illumination (photons/s/pixel)
  m = re.search(r'^ILLUMINATION:\s*(\S+)', line)
  if m: I = float(m.group(1))
  # QE (Illumination * QE = current, e/s/pixel)
  m = re.search(r'^QE:\s*(\S+)', line)
  if m: QE = float(m.group(1))

  # RNG seed
  m = re.search(r'^RNGSEED:\s*(\d+)', line)
  if m: rngseed = int(m.group(1))

  # Noise
  m = re.search(r'^NOISE:\s*(\S+)\s+(\S+)', line)
  if m:
    noisemode = m.group(1)
    if noisemode != 'none':
      noisefile = m.group(2)

  # Wavelength mode, affects BFE choice
  m = re.search(r'^WAVEMODE:\s*(\S+)', line)
  if m: wavemode = m.group(1)

  # BFE
  m = re.search(r'^BFE:\s*(\S+)', line)
  if m: bfemode = m.group(1)

  # linear IPC
  m = re.search(r'^L_IPC:\s*(\S+)\s+(\S.*)', line)
  if m:
    lipcmode = m.group(1)
    if lipcmode == 'true':
      lipc_alpha_str = m.group(2).split(" ")
      lipc_alpha = [ float(lipc_alpha_str[x]) for x in range(len(lipc_alpha_str)) ]

  # non-linearity beta
  m = re.search(r'^NL:\s*(\S+)\s+(\S.*)', line)
  if m:
    nlmode = m.group(1)
    if nlmode == 'quadratic':
      nlbeta = float(m.group(2))
    elif nlmode == 'quartic':
      nlcoeffs_str = m.group(2).split(" ")
      nlcoeffs_arr = [ float(nlcoeffs_str[x]) for x in range(len(nlcoeffs_str)) ]
  # Reset level (in e)
  m = re.search(r'^RESET_E:\s*(\S+)', line)
  if m: resetlevel = float(m.group(1))

  # Output file
  m = re.search(r'^OUTPUT:\s*(\S+)', line)
  if m: outfile = m.group(1)

# data cube attributes
N = nx = ny = get_nside(formatpars) # possibly redundant with nx,ny
# Reference pixels hard-coded to 4 rows/cols around border, true for
# all except WFC3 which has 5
xmin,xmax,ymin,ymax = 4,N-4,4,N-4 # Extent of non-reference pixels
nt_step = tsamp*substep # number of tot timesteps depending on convergence needs
delta_t = (delta_tsamp*tsamp)/nt_step # time between timesteps
allQ = np.zeros((substep, nx, ny))
data_cube_Q = np.zeros((tsamp, nx, ny))
data_cube_S = np.zeros_like(data_cube_Q)
offset_frame = np.zeros((1, nx, ny))
count = 1
reset_count = 0

print('side length =',N)
print('samples:', tsamp, 'x', delta_tsamp, 's; # substep =', substep)
print('Illumination:', I, 'ph/s/pix; QE =', QE)

print('RNG seed ->', rngseed)
numpy.random.seed(rngseed)

# Define function to build distribution of pairs of points
# myMean is mean number of *photons* per pixel (number of e is 2x larger)
def PairPoisson(myMean, myShape):
  Ptot = np.zeros(myShape)
  (yg,xg) = myShape
  for j in range(2*QY_offset+1):
    for i in range(2*QY_offset+1):
      dP = np.random.poisson(myMean*QY_p2[j,i], (yg+2*QY_offset, xg+2*QY_offset))
      Ptot += dP[QY_offset:yg+QY_offset, QY_offset:xg+QY_offset]
      Ptot += dP[j:yg+j, i:xg+i]
  return(Ptot)

# Reset first frame if needed
offset_frame[:,:,:] = resetlevel
if 0 in reset_frames:
  allQ[0,:,:] = data_cube_Q[0,:,:] = offset_frame

# Start with 0 charge in the first frame (t=0)
mean = I*delta_t

for tdx in range(1, nt_step):
  # Create realization of charge
  # This version uses less memory, but probably still sub-optimal
  idx = tdx%substep
  # Charge accumulates, dependent on quantum efficiency of the pixel
  # and the brighter-fatter effect.  First timestep is Poisson realization
  if (tdx==1):
    allQ[idx,:,:] = allQ[idx-1,:,:]
    allQ[idx,xmin:xmax,ymin:ymax] += np.random.poisson(
      QE*mean*(1.-QY_omega)/(1.+QY_omega), allQ[idx,xmin:xmax,xmin:xmax].shape)
    allQ[idx,xmin:xmax,ymin:ymax] += PairPoisson(QE*mean*QY_omega/(1.+QY_omega), allQ[idx,xmin:xmax,xmin:xmax].shape)
  else:
    # If not the first step, and the brighter-fatter effect is turned
    # on, then now loop through all pixels
    if bfemode=='true':
      # Calculate the area defect by taking a convolution of the bfe
      # kernel (flipped in the calc_area_defect function) and the charge Q
      # Area defects are magnified by (1+QY_omega)/(1-QY_omega) since I am only applying them
      # to the 1-electron events
      if wavemode=='ir':
        a_coeff = get_bfe_kernel_5x5_ir()
      elif wavemode=='vis':
        a_coeff = get_bfe_kernel_5x5_vis()
      else:
        print("wavemode set to unknown value, defaulting BFE to IR")
        a_coeff = get_bfe_kernel_5x5()
      area_defect = calc_area_defect(
        (1.+QY_omega)/(1.-QY_omega)*a_coeff, allQ[idx-1,xmin:xmax,ymin:ymax])
      meanQ = area_defect*mean*QE*(1.-QY_omega)/(1.+QY_omega)
      allQ[idx,xmin:xmax,ymin:ymax] = allQ[idx-1,xmin:xmax,ymin:ymax] + \
          np.random.poisson(meanQ)
      allQ[idx,xmin:xmax,ymin:ymax] += PairPoisson(QE*mean*QY_omega/(1.+QY_omega), allQ[idx,xmin:xmax,xmin:xmax].shape)
    else:
      # Otherwise Poisson draw the charge as before
      allQ[idx,:,:] = allQ[idx-1,:,:]
      allQ[idx,xmin:xmax,ymin:ymax] += np.random.poisson(
        QE*mean*(1.-QY_omega)/(1.+QY_omega), allQ[idx,xmin:xmax,ymin:ymax].shape)
      allQ[idx,xmin:xmax,ymin:ymax] += PairPoisson(QE*mean*QY_omega/(1.+QY_omega), allQ[idx,xmin:xmax,xmin:xmax].shape)
  if (idx==0):
    data_cube_Q[count,:,:] = allQ[idx,:,:]
    allQ = np.zeros((substep, nx, ny))
    # if this is a reset frame set start to offset
    if count in reset_frames:
      allQ[0,:,:] = offset_frame
    else:
      allQ[0,:,:] = data_cube_Q[count,:,:]

    count += 1
    print("time: %d" %count)
    
# Add in IPC before the noise if the mode is turned on
if (lipcmode=='true'):
  data_cube_Q[:,xmin:xmax,ymin:ymax] = calculate_ipc(
    data_cube_Q[:,xmin:xmax,ymin:ymax], lipc_alpha)
else:
  pass

# Apply non-linearity if mode turned on; assumed to act after IPC
if (nlmode=='quadratic'):
  data_cube_Q[:,xmin:xmax,ymin:ymax] -= (1.E-6*nlbeta) * \
      data_cube_Q[:,xmin:xmax,ymin:ymax]**2
  print("Applying non-linearity at leading order coefficient (quadratic term)")
elif (nlmode=='quartic'):
  data_cube_Q[:,xmin:xmax,ymin:ymax] += 1.E-6*nlcoeffs_arr[0] * \
      data_cube_Q[:,xmin:xmax,ymin:ymax]**2 + 1.E-6*nlcoeffs_arr[1] * \
      data_cube_Q[:,xmin:xmax,ymin:ymax]**3 + 1.E-6*nlcoeffs_arr[2] * \
      data_cube_Q[:,xmin:xmax,ymin:ymax]**4
  print("Applying non-linearity polynomial to quartic term")
else:
  print("No additional non-linearity (Beta) applied")
  pass

# Read in the read noise from a fits file generated with Bernie's ngxhrg
# noisemode 'last' uses one realization because the full one takes
# a long time to create
if noisemode == 'last':
  noise = fitsio.read(noisefile)
  data_cube_Q[-1,:,:] += noise  # Adding to only the final time
elif noisemode == 'full':
  noise = fitsio.read(noisefile)
  data_cube_Q += noise  # Adding the noise at all reads
else:
  # a small amount of Gaussian noise, 12 e
  data_cube_Q += 10*normal(size=(tsamp,N,N))

# Convert charge to signal, clipping values<0 and >2**16
data_cube_S = np.array(
  np.clip(data_cube_Q/gain, 0, 65535), dtype=np.uint16)

# Write simple header, todo: add more thorough comments
hdr = FITSHDR()
hdr['GAIN'] = gain
hdr['ILLUMIN'] = I
hdr['QE'] = QE
hdr['RNGSEED'] = rngseed
if (lipcmode=='true'):
  hdr['LINIPC'] = lipc_alpha[0]
if (nlmode=='quadratic'):
  hdr['BETA'] = nlbeta
if (nlmode=='quartic'):
  hdr['NLCOEFFS_c2'] = nlcoeffs_arr[0]
  hdr['NLCOEFFS_c3'] = nlcoeffs_arr[1]
  hdr['NLCOEFFS_c4'] = nlcoeffs_arr[2]
if (bfemode=='true'):
  hdr['BFE_A00'] = a_coeff[2][2]  # Hard-coded to expect 5x5 a coeffs

# Open up an example DCL flat file and save the data cube
#dclfile = 'Set_001_Test_0002.fits'
fitsio.write(outfile, data_cube_S, header=hdr, clobber=True)

# Try compression of data cube into file
# End of script
