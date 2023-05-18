import numpy
from astropy.io import fits
import fitsio

### INPUT/OUTPUT ###

# read a segment of a file
def read_segment(infile, bounds=None):
  use_fitsio = True
  xmin = 0; xmax = 4096; ymin = 0; ymax = 4096
  if bounds!=None: [ymin,ymax,xmin,xmax] = bounds
  map = numpy.zeros((ymax-ymin, xmax-xmin))
  if use_fitsio:
    fileh = fitsio.FITS(infile)
    if len(fileh)==1:
      map[:,:] = fileh[0][ymin:ymax,xmin:xmax]
    if len(fileh)==2:
      nt = int((fileh[1].read_header())['NAXIS3'])
      map[:,:] = fileh[1][0,0,ymin:ymax,xmin:xmax].astype(numpy.float64) - fileh[1][0,nt-1,ymin:ymax,xmin:xmax].astype(numpy.float64)
    fileh.close()
  else:
    hdul = fits.open(infile)
    if len(hdul)==1:
      map[:,:] = hdul[0].data[ymin:ymax,xmin:xmax]
    if len(hdul)==2:
      map[:,:] = hdul[1].data[0,0,ymin:ymax,xmin:xmax].astype(numpy.float64) - hdul[1].data[0,-1,ymin:ymax,xmin:xmax].astype(numpy.float64)
    hdul.close()
  return(map)

### COVARIANCE MATRIX AND TEMPLATE MANAGEMENT ###

nTemplate = 6

# builds a covariance matrix C for an (ny,nx)-sized postage stamp
#
# inputs:
#  size -> tuple (ny,nx) size of array
#  powerspec -> 2D array of the power spectrum of size (nyPS, nxPS), both even
#    powerspec[jy,jx] corresponds to jy/nyPS, jx/nxPS cycles per pixel
#    (will crash if nyPS,nxPS <= 2 ny,2 nx -- results in pattern wrapping)
#    normalization: sum is the total variance
#  branch -> tuple (alpha,beta) of integers, describing offsets:
#    negative integer for data with no offset (default)
#
# output covariance matrix:
#   output[j1,j2] = d(data_j1)/d(xi_j1,alpha), d(data_j2)/d(xi_j2,beta)
#   (shaped as ny*nx x ny*nx matrix)
# 
def getCovar(size, powerspec, branch=(-1,-1)):
  (ny,nx)=size
  (nyPS,nxPS) = numpy.shape(powerspec)
  ps = numpy.copy(powerspec.astype(numpy.complex128)) # convert to complex

  # checks
  if nyPS<=2*ny or nxPS<=2*nx:
    print('astromutils.getCovar: size error: asked for subarray of size', size, 'from gridded power spectrum of size', numpy.shape(powerspec))
    exit()

  kx = numpy.zeros((nyPS,nxPS))
  for i in range(-(nxPS//2),nxPS//2):
    kx[:,i] = 2*numpy.pi*i/float(nxPS)
  ky = numpy.zeros((nyPS,nxPS))
  for i in range(-(nyPS//2),nyPS//2):
    ky[i,:] = 2*numpy.pi*i/float(nyPS)

  # counterintutive that we want to put the +ik on the alpha side and -ik on the beta
  # it is because the formula is
  # Cov[I(x'[alpha]), I(x[beta])] = sum P(k) e^{ikx} e^{-ikx'} Talpha*(k) Tbeta(k)
  # *but* fft2 is 'forward' and introduces another complex conjugate.
  # (ifft2 has the wrong normalization by default)

  # note template 0 doesn't do anything (repeat of -1)

  # alpha side
  if branch[0]==1: ps = 1j*ps*kx
  if branch[0]==2: ps = 1j*ps*ky
  if branch[0]==3: ps = -.5*ps*kx**2
  if branch[0]==4: ps = -1*ps*kx*ky
  if branch[0]==5: ps = -.5*ps*ky**2

  # beta side
  if branch[1]==1: ps = -1j*ps*kx
  if branch[1]==2: ps = -1j*ps*ky
  if branch[1]==3: ps = -.5*ps*kx**2
  if branch[1]==4: ps = -1*ps*kx*ky
  if branch[1]==5: ps = -.5*ps*ky**2

  # build the correlation function
  CorrFunc = numpy.fft.fftshift(numpy.real(numpy.fft.fft2(ps))).astype(numpy.float64)
  # central entry is in CorrFunc[nyPS//2,nxPS//2]

  # now make covariance out of this
  cov = numpy.zeros((ny,nx,ny,nx))
  for j in range(ny):
    for i in range(nx):
      cov[j,i,:,:] = CorrFunc[nyPS//2-j:nyPS//2-j+ny,nxPS//2-i:nxPS//2-i+nx]
  return(cov.reshape((ny*nx,ny*nx)))

# fast symmetric multiplication and C-inverse for big matrices

def force_sym_mult(A,B):
  (N,N1) = numpy.shape(A)
  m = 4
  C = numpy.zeros((N,N))
  for i in range(m):
    imin = (i*N)//m; imax = ((i+1)*N)//m
    for j in range(i+1):
      jmin = (j*N)//m; jmax = ((j+1)*N)//m
      C[imin:imax,jmin:jmax] = A[imin:imax,:]@B[:,jmin:jmax]
      if j<i: C[jmin:jmax,imin:imax] = C[imin:imax,jmin:jmax].T
  return(C)

def syminv(A):
  (N,N1) = numpy.shape(A)
  N2 = N//2
  if N1!=N:
    print('Error in syminv -> wrong shape.')
    exit()

  # Strassen algorithm
  R1 = numpy.linalg.inv(A[:N2,:N2])
  R2 = A[N2:,:N2]@R1
  R3 = R2.T
  R6 = numpy.linalg.inv(A[N2:,:N2]@R3-A[N2:,N2:])

  # output
  C = numpy.zeros_like(A)
  C[N2:,:N2] = R6@R2
  C[:N2,N2:] = C[N2:,:N2].T
  C[:N2,:N2] = R1-force_sym_mult(R3,C[N2:,:N2])
  C[N2:,N2:] = -R6
  return(C) 

