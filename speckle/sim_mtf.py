import sys
import numpy
import re
from astropy.io import fits
import fitsio
from fitsio import FITS,FITSHDR

nside = 4096 # nside
pad = 128
N = nside+pad

Amp = 0.
gain = [1.]
dir = ''
slitdim = [0.,0.,0.]
rngseed = 0
cd = 0.
alpha = 0.
cdsnoise = 18.
outfile = ''
nrepeat = 1
custom = 0
offsetPattern = 0

nPerturbations = 6

# read configuration file
with open(sys.argv[1]) as f: config = f.readlines()
for l in config:

  # number of realizations
  m = re.search(r'^REPEAT\:\s*(\S+)', l)
  if m: nrepeat = int(m.group(1))

  # amplitude
  m = re.search(r'^AMPLITUDE\:\s*(\S+)', l)
  if m: Amp = float(m.group(1))

  # seed
  m = re.search(r'^RNGSEED\:\s*(\S+)', l)
  if m: rngseed = int(m.group(1))

  # output
  m = re.search(r'^OUTFILE\:\s*(\S+)', l)
  if m: outfile = m.group(1)

  # gain
  m = re.search(r'^GAIN\:\s*(\S.*)', l)
  if m:
    gain = re.split(r'\s+', m.group(1))
    for j in range(len(gain)): gain[j] = float(gain[j])

  # IPC
  m = re.search(r'^ALPHA\:\s*(\S+)', l)
  if m: alpha = float(m.group(1))

  # direction, H or V (lower case OK)
  m = re.search(r'^DIR\:\s*(\w)', l)
  if m: dir = m.group(1).upper()

  # slit dimensions, projected onto detector
  # first the long dimension, then the short dimension, then the center spacing
  m = re.search(r'^SLITDIM\:\s*(\S.*)', l)
  if m:
    slitdim = re.split(r'\s+', m.group(1))
    for j in range(len(slitdim)): slitdim[j] = float(slitdim[j])

  # custom information
  m = re.search(r'^CUSTOM\:\s*(\S.*)', l)
  if m:
    customSettings = re.split(r'\s+', m.group(1))
    custom = int(customSettings[0])
    for k in range(1,len(customSettings)): customSettings[k] = float(customSettings[k])

  # custom information
  m = re.search(r'^OFFSET\:\s*(\S.*)', l)
  if m:
    offsetSettings = re.split(r'\s+', m.group(1))
    offsetPattern = int(offsetSettings[0])
    for k in range(1,len(offsetSettings)): offsetSettings[k] = float(offsetSettings[k])

# realization
numpy.random.seed(rngseed)

# Checking, settings
nc = len(gain) # number of column groups of different gain
if custom>0:
  print('Customized case '+str(custom)+', parameters:', customSettings[1:])
  dir = 'X'; slitdim = [300,100,100] # avoid error for custom slit
if dir=='':
  print('Error: no DIR specified.'); exit()
if len(slitdim)!=3:
  print('Error: wrong SLITDIM length:', len(slitdim), ' <-- should be 3'); exit()
if rngseed==0:
  print('Error: no RNG seed specified.'); exit()
if outfile=='':
  print('Error: no output specified.'); exit()

#
# oversampling
cmax = 2*numpy.amax([slitdim[0], slitdim[1]+slitdim[2]])
ov = 1
while (cmax>=N*ov): ov*=2

print('GENERAL SETTINGS:\n')
print('Amplitude:', Amp, 'e')
print('Gain:', gain, 'e/DN')
print('IPC alpha:', alpha)
print('Direction:', dir)
print('Slit dimensions:', slitdim, 'cy/array')
print('using oversampling:', ov, 'x')
print('random seed:', rngseed)
print('output:', outfile, '->', nrepeat, 'realizations')
print('')

# perturbations
imOffset = numpy.zeros((nPerturbations,nside,nside))
kx = numpy.zeros((N*ov,N*ov))
ky = numpy.zeros((N*ov,N*ov))
for i in range(-(N*ov//2),N*ov//2):
  kx[:,i] = 2*numpy.pi*i/float(N*ov)
  ky[i,:] = 2*numpy.pi*i/float(N*ov)

if offsetPattern==1:
  for i in range(512):
    for j in range(512):
      imOffset[1,8*j:8*(j+1),8*i:8*(i+1)] = -.02
      if (i+j)%2==1: imOffset[0,8*j:8*(j+1),8*i:8*(i+1)] = .02
if offsetPattern==2:
  for i in range(512):
    for j in range(512):
      imOffset[1,8*j:8*(j+1),8*i:8*(i+1)] = -.02
      if (i+j)%2==1: imOffset[1,8*j:8*(j+1),8*i:8*(i+1)] = .02
  for i in range(32768):
    xc = numpy.random.randint(4,4096); yc = numpy.random.randint(4,4096)
    kmax = numpy.random.randint(1,16)
    for k in range(kmax):
      imOffset[2,yc-1,xc-k] = -.02
      imOffset[2,yc,xc-k] = .02
    xc = numpy.random.randint(4,4096); yc = numpy.random.randint(4,4096)
    kmax = numpy.random.randint(1,16)
    for k in range(kmax):
      imOffset[2,yc-k,xc-1] = -.02
      imOffset[2,yc-k,xc] = .02
  imOffset[3,:,:] = imOffset[1,:,:]**2
  imOffset[5,:,:] = imOffset[2,:,:]**2
  for i in range(65536):
    xc = numpy.random.randint(8,4088); yc = numpy.random.randint(8,4088)
    imOffset[0,yc-1:yc+2,xc-1:xc+2] = -.02
  for a in range(1,6):
    imOffset[a,:,:] *= 1.+imOffset[0,:,:]
imOffset[:,:,:4] = 0.; imOffset[:,:,-4:] = 0.; imOffset[:,:4,:] = 0.; imOffset[:,-4:,:] = 0

# loop over realizations
for krepeat in range(1,nrepeat+1):
  print('generating '+outfile+'_{:03d}.fits'.format(krepeat),'...')

  # make E-field map
  E_FT = numpy.zeros((ov*N,ov*N), dtype=numpy.complex128)

  if custom==0:
    # double slit
    L = int(slitdim[0]*N/nside)
    W = int(slitdim[1]*N/nside)
    sp = int(slitdim[2]*N/nside)
    E_FT[:W,:L] = numpy.random.normal(size=(W,L)) + 1j*numpy.random.normal(size=(W,L))
    E_FT[sp:sp+W,:L] = numpy.random.normal(size=(W,L)) + 1j*numpy.random.normal(size=(W,L))
  if custom==1:
    # triangle array
    # customSettings = [1, width, trianglesidelength]
    L = int(customSettings[1]*N/nside)
    for j in range(3):
      theta = 2.*numpy.pi * (j/3. + 1./24.)
      xc = int(customSettings[2]*N/nside*(1+numpy.cos(theta))/numpy.sqrt(3))
      yc = int(customSettings[2]*N/nside*(1+numpy.sin(theta))/numpy.sqrt(3))
      E_FT[yc:yc+L,xc:xc+L] = numpy.random.normal(size=(L,L)) + 1j*numpy.random.normal(size=(L,L))
      for jj in range(L):
        for ii in range(L):
          if (jj-L//2)**2+(ii-L//2)**2>=(L//2)**2:
            E_FT[yc+jj,xc+ii] = 0.
  if custom==2:
    # tilted slit
    # customSettings = [2, width, length, angle]
    W = int(customSettings[1]*N/nside)
    L = int(customSettings[2]*N/nside)
    T = 2*(L+W)
    E_FT[:T,:T] = numpy.random.normal(size=(T,T)) + 1j*numpy.random.normal(size=(T,T))
    ca = numpy.cos(customSettings[3]*numpy.pi/180)
    sa = numpy.sin(customSettings[3]*numpy.pi/180)
    for j in range(T):
      for i in range(T):
        xx =  (i-L-W)*ca + (j-L-W)*sa
        yy = -(i-L-W)*sa + (j-L-W)*ca
        if (xx<-L) or (xx>L) or (yy<-W) or (yy>W): E_FT[j,i] = 0.
  if custom==3:
    # hexagon array
    # customSettings = [1, width, trianglesidelength]
    L = int(customSettings[1]*N/nside)
    for j in range(7):
      theta = 2.*numpy.pi * (j/6. + 1./24.)
      xc = int(customSettings[2]*N/nside*(1+numpy.cos(theta))/numpy.sqrt(3))
      yc = int(customSettings[2]*N/nside*(1+numpy.sin(theta))/numpy.sqrt(3))
      if j==6: xc = yc = int(customSettings[2]*N/nside)
      E_FT[yc:yc+L,xc:xc+L] = numpy.random.normal(size=(L,L)) + 1j*numpy.random.normal(size=(L,L))

  # auto-convolve
  E = numpy.fft.ifft2(E_FT)
  Intensity = numpy.fft.fft2(numpy.abs(E)**2)

  #
  # smooth with pixel tophat & charge diffusion
  SmoothVec = numpy.zeros((ov*N,))
  for i in range(-ov*(N//2)+1,ov*N//2):
    kpix = 2.*numpy.pi*float(i)/N
    if i==0:
      SmoothVec[i] = 1.
    else:
      SmoothVec[i] = numpy.sin(kpix/2.)/(kpix/2.)
    SmoothVec[i] *= numpy.exp(-cd**2*kpix**2/2)
  for i in range(ov*N):
    Intensity[:,i] *= SmoothVec[i]
    Intensity[i,:] *= SmoothVec[i]

  # take out a part of the image
  Intensity0 = numpy.copy(Intensity)
  Intensity = numpy.real(numpy.fft.ifft2(Intensity))
  C = Amp/numpy.mean(Intensity)
  Intensity *= C
  Intensity0 *= C
  if dir=='H':
    Intensity = numpy.transpose(Intensity)
    Intensity0 = numpy.transpose(Intensity0)

  # down-sample and clip
  ImageNoNoise = Intensity[0:ov*nside:ov,0:ov*nside:ov]
  # include image distortions
  ImageNoNoise = ImageNoNoise*(1.+imOffset[0,:,:]) + (
                    imOffset[1,:,:]*numpy.real(numpy.fft.ifft2(1j*Intensity0*kx))[0:ov*nside:ov,0:ov*nside:ov]
                   +imOffset[2,:,:]*numpy.real(numpy.fft.ifft2(1j*Intensity0*ky))[0:ov*nside:ov,0:ov*nside:ov]
                   +imOffset[3,:,:]*numpy.real(numpy.fft.ifft2(-.5*Intensity0*kx*kx))[0:ov*nside:ov,0:ov*nside:ov]
                   +imOffset[4,:,:]*numpy.real(numpy.fft.ifft2(-1*Intensity0*kx*ky))[0:ov*nside:ov,0:ov*nside:ov]
                   +imOffset[5,:,:]*numpy.real(numpy.fft.ifft2(-.5*Intensity0*ky*ky))[0:ov*nside:ov,0:ov*nside:ov]
                 )
  ImageNoNoise = numpy.where(ImageNoNoise>0,ImageNoNoise,0.)

  # add noise, IPC
  ImagePoisson = numpy.random.poisson(ImageNoNoise).astype(numpy.float64)
  ImageIPC = numpy.copy(ImagePoisson)
  ImageIPC[5:-4,4:-4] += alpha*(ImagePoisson[4:-5,4:-4]-ImagePoisson[5:-4,4:-4])
  ImageIPC[4:-5,4:-4] += alpha*(ImagePoisson[5:-4,4:-4]-ImagePoisson[4:-5,4:-4])
  ImageIPC[4:-4,5:-4] += alpha*(ImagePoisson[4:-4,4:-5]-ImagePoisson[4:-4,5:-4])
  ImageIPC[4:-4,4:-5] += alpha*(ImagePoisson[4:-4,5:-4]-ImagePoisson[4:-4,4:-5])
  #
  # reference pixels
  ImageIPC[:4,:] = 0.
  ImageIPC[-4:,:] = 0.
  ImageIPC[:,:4] = 0.
  ImageIPC[:,-4:] = 0.
  ImageNoise = ImageIPC + numpy.random.normal(size=(nside,nside))*numpy.sqrt(cdsnoise)

  # gain
  for i in range(nside): ImageNoise[:,i] /= gain[(i*nc)//nside]

  # quantize
  ImageNoise = numpy.where(ImageNoise>0,ImageNoise,0.)
  ImageNoise = numpy.where(ImageNoise<65535,ImageNoise,65535.)
  ImageOut = numpy.rint(ImageNoise).astype(numpy.int32)

  # output
  hdr = FITSHDR()
  hdr['GAIN'] = gain
  hdr['RNGSEED'] = rngseed
  hdr['IPCALPHA'] = alpha
  fitsio.write(outfile+'_{:03d}.fits'.format(krepeat), ImageOut, header=hdr, clobber=True)

# offset "key"
hdr = FITSHDR()
fitsio.write(outfile+'_offsetkey.fits', imOffset.astype(numpy.float32), header=hdr, clobber=True)
