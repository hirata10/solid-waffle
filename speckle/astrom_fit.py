import sys
import time
import numpy
import re
from astropy.io import fits
from scipy.ndimage.filters import convolve
from scipy import ndimage
import astromutils

startInput = 1
inputdata = []
nrepeat = []
dir = []
cppinit = []

outstem = ''

nside = 4096

# (20,.25)

niter = 20 # number of iterations for offset map
condition_factor = .5

QuadIter=True
ChebIter=False

if QuadIter:
  lmax = 6.
  condition_factor = 2.6833/lmax
  condition_gamma = 1.8334/lmax

if ChebIter:
  # reduce by factor of Tn(1/alphaCheb)
  lmax = 6.
  nCheb = 3
  kCheb = 1 # relatively prime to nCheb
  alphaCheb = .75
  niter = (niter//nCheb)*nCheb

ForceFlat=True
flatprior = 1.0e-4
bgSmooth = 3

cflat = None # no input flat field
use_subregion = False

# only compute diagonal blocks in Fisher
DiagFisher = True

buildRange = None

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

  # File number to start
  m = re.search(r'^START:\s*(\d+)', line)
  if m: startInput = int(m.group(1))

  # Bin sizes
  m = re.search(r'^NBIN:\s*(\d+)\s+(\d+)', line)
  if m:
    nx = int(m.group(1))
    ny = int(m.group(2))

  # Which to build?
  m = re.search(r'^BUILD:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
  if m: buildRange = [int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))]

  m = re.search(r'^OUTPUT:\s*(\S+)', line)
  if m: outstem = m.group(1)

  m = re.search(r'^CFLAT:\s*(\S+)', line)
  if m: cflat = m.group(1)

  m = re.search(r'^POWERSUB', line)
  if m: use_subregion = True

# end read config

# get a particular row
if len(sys.argv)>2:
  thisrow = int(sys.argv[2])
  outstem += '_Row{:02d}'.format(thisrow)
  buildRange[0] = thisrow
  buildRange[1] = thisrow+1
  if len(sys.argv)>3: buildRange[1] = int(sys.argv[3]) + buildRange[0]

if outstem=='':
  print('Error: no output specified.')
  exit()

n_input_group = len(inputdata)
max_nrepeat = numpy.amax(numpy.asarray(nrepeat))
dx = nside//nx; dy = nside//ny
SCAMask = numpy.zeros((n_input_group, nside, nside), dtype=numpy.uint16) # 16 bit mask
SCAMask[:,:,:4] = 1; SCAMask[:,:,-4:] = 1; SCAMask[:,:4,:] = 1; SCAMask[:,-4:,:] = 1 # reference pixels

# Tell the user what we learned
print('Using files:')
for k in range(n_input_group): print('   ',inputdata[k], 'repeats:', nrepeat[k])
print('max. number of repeats =', max_nrepeat)
print('number of superpixels: {:d} (in y) x {:d} (in x) = {:d} (total)'.format(ny,nx,nx*ny))
print('flat fields in', cflat)

# Check build ranges
if buildRange is None: buildRange = [0,ny,0,nx]
print('building range: {:3d}<=ty<{:3d}, {:3d}<=tx<{:3d}'.format(buildRange[0],buildRange[1],buildRange[2],buildRange[3]))
if buildRange[0]<0 or buildRange[1]>ny or buildRange[2]<0 or buildRange[3]>nx:
  print('Error: illegal build range.')
  exit()

# get the cflat-field
cflat_image = numpy.ones((4096,4096))
if cflat:
  with fits.open(cflat) as hdul: cflat_image = hdul[0].data
  cflat_image2 = numpy.copy(cflat_image)
  cflat_image2[:4,:] = 1.
  cflat_image2[-4:,:] = 1.
  cflat_image2[:,:4] = 1.
  cflat_image2[:,-4:] = 1.
  for k in range(n_input_group):
    SCAMask[k][numpy.where(cflat_image2<0.5)] |= 2
    SCAMask[k][numpy.where(cflat_image2>1.5)] |= 2
    # now mask 8-pixel region around the defect
    block = numpy.where(numpy.bitwise_and(SCAMask[k,4:-4,4:-4],2)>0, 1, 0).astype(numpy.int16)
    block2 = numpy.zeros((4096,4096), dtype=numpy.int16)
    for ddy in [-1,0,1]:
      for ddx in [-1,0,1]:
        block2[4+ddy:-4+ddy,4+ddx:-4+ddx] += block
    SCAMask[k][numpy.where(block2>0)] |= 4
    del block; del block2
  cflat_image = numpy.where(numpy.bitwise_or(cflat_image<0.5,cflat_image>1.5), 1., cflat_image)

# Tell us about the mask
print('bit true count in mask:')
for x in range(3): print('    bit', x, numpy.sum(1-numpy.prod(numpy.where(numpy.bitwise_and(SCAMask,1<<x)==0,1,0),axis=0)))
print('')

# Get the power spectrum of each class of files
downsample = 4 # down-sampling factor of power spectrum
PS = numpy.zeros((n_input_group,nside,nside))
PS_downsample = numpy.zeros((n_input_group,downsample,downsample,nside//downsample,nside//downsample))
norm = numpy.zeros((n_input_group, max_nrepeat))
ntot = numpy.sum(nrepeat)
j=0
MM = 4
qnorm = numpy.zeros((ntot,MM*ny,MM*nx))
for k in range(n_input_group):
  for ik in range(nrepeat[k]):
    print('FT FILE:', re.sub(r'\*', '{:03d}'.format(ik+startInput), inputdata[k]))
    map = astromutils.read_segment(re.sub(r'\*', '{:03d}'.format(ik+startInput), inputdata[k]))/cflat_image
    # normalize, remove zero mode & reference pixels
    norm[k,ik] = numpy.mean(map[4:-4,4:-4])
    map[4:-4,4:-4] /= norm[k,ik]
    map[4:-4,4:-4] -= 1.
    map[:4,:] = 0.; map[-4:,:] = 0.; map[:,:4] = 0.; map[:,-4:] = 0.
    PS[k,:,:] = PS[k,:,:] + numpy.abs(numpy.fft.fft2(map)**2)/nrepeat[k]

    for yblock in range(downsample):
      for xblock in range(downsample):
        PS_downsample[k,yblock,xblock,:,:] = PS_downsample[k,yblock,xblock,:,:] +\
          numpy.abs(numpy.fft.fft2(map[yblock*nside//downsample:(yblock+1)*nside//downsample,\
          xblock*nside//downsample:(xblock+1)*nside//downsample])**2)/nrepeat[k]

    for ty in range(MM*ny):
      ymin = ty*dy//MM; ymax = ymin+dy//MM
      for tx in range(MM*nx):
        xmin = tx*dx//MM; xmax = xmin+dx//MM
        qnorm[j,ty,tx] = numpy.mean(map[ymin:ymax,xmin:xmax])
    j+=1
del map
PS_downsample *= (downsample/float(nside))**4 # new version of power spectrum

qnorm = numpy.mean(qnorm, axis=0)
qnorm = numpy.median(qnorm.reshape((ny,MM,nx,MM)), axis=(1,3))
qnorm = ndimage.median_filter(qnorm, size=(bgSmooth,bgSmooth), mode='reflect')
print('normalization factors -->', numpy.amin(qnorm), numpy.median(qnorm), numpy.amax(qnorm))
#
# write normalization map to a file
hdu = fits.PrimaryHDU(qnorm.astype(numpy.float32))
hdu.writeto(outstem+'_qnorm.fits', overwrite=True)

# Write the power spectrum to a file
hdu = fits.PrimaryHDU(numpy.fft.fftshift(PS, axes=(1,2)).astype(numpy.float32))
hdu.writeto(outstem+'_powerspec.fits', overwrite=True)
#
# smeared version
smooth_length = 9 # make odd
PS_smooth = numpy.zeros((n_input_group,nside//downsample,nside//downsample))
for k in range(n_input_group):
  PS_smooth[k,:,:] = convolve(PS[k,:,:],numpy.full((smooth_length,smooth_length),1./smooth_length**2),mode='wrap')[::downsample,::downsample]
PS_smooth *= downsample**2/nside**4
hdu = fits.PrimaryHDU(numpy.log(numpy.fft.fftshift(PS_smooth, axes=(1,2))).astype(numpy.float32))
hdu.writeto(outstem+'_lnPk_smooth.fits', overwrite=True)
for k in range(n_input_group):
  var = numpy.sum(PS_smooth[k,:,:])
  cy = ny//downsample
  cx = nx//downsample
  varsub = numpy.sum(PS_smooth[k,cy:-cy,cx:-cx])
  print('total variance in group {:2d} is {:9.6f}, per superpix {:9.6f}'.format(k,var,varsub))

# and dx map
dxmap = numpy.zeros((astromutils.nTemplate, nside, nside))
# and dx map with filtering
sigmaprior = numpy.asarray([1e-4,.1,.1,1e-4,1e-4,1e-4])
dxfiltmap = numpy.zeros((astromutils.nTemplate, nside, nside))

# mask here
hdu = fits.PrimaryHDU(SCAMask)
hdu.writeto(outstem+'_mask.fits', overwrite=True)

# change to ty in range(ny), tx in range(nx) for full coverage
for ty in range(buildRange[0], buildRange[1]):
  ymin = ty*dy; ymax = ymin+dy
  for tx in range(buildRange[2], buildRange[3]):
    xmin = tx*dx; xmax = xmin+dx

    # get the weighting covariance matrices
    yblock = ty*dy//(nside//downsample)
    xblock = tx*dx//(nside//downsample)
    print('super-pixel', ty, tx, yblock, xblock)
    if use_subregion: PS_smooth = PS_downsample[:,yblock,xblock,:,:]
    CovWT = numpy.zeros((n_input_group,dy*dx,dy*dx))
    for k in range(n_input_group):
      print('Computing covariance matrix for group {:2d} ...'.format(k))
      CovWT[k,:,:] = astromutils.getCovar((dy,dx), PS_smooth[k,:,:])

    for cy in range(niter):
      Fisher = numpy.zeros((astromutils.nTemplate*dy*dx,astromutils.nTemplate*dy*dx))
      any_good = numpy.zeros((dy*dx,), dtype = numpy.bool_)

      dNLLdp = numpy.zeros((astromutils.nTemplate, dy, dx))
      for a in range(astromutils.nTemplate): dNLLdp[a,:,:] = -dxmap[a,ymin:ymax,xmin:xmax]/sigmaprior[a]**2

      for k in range(n_input_group):
        print('group {:2d} [{:4d}:{:4d},{:4d}:{:4d}] cy{:4d} at {:s}'.format(k,ymin,ymax,xmin,xmax,cy,time.asctime(time.localtime(time.time()))))

        # C,a-templates for this group
        dCdp = numpy.zeros((astromutils.nTemplate, dy*dx, dy*dx))
        for a in range(astromutils.nTemplate): dCdp[a,:,:] = astromutils.getCovar((dy,dx), PS_smooth[k,:,:], branch=(a,-1))
        # (will add more later)
        # get the C^-1
        thisC = numpy.copy(CovWT[k,:,:])
        if cy>0:
          xi = dxmap[:,ymin:ymax,xmin:xmax].reshape((astromutils.nTemplate,dy*dx))
          for a in range(astromutils.nTemplate):
            d_thisC = astromutils.getCovar((dy,dx), PS_smooth[k,:,:], branch=(a,-1)) * numpy.tile(xi[a,:],(dy*dx,1)).T
            thisC += d_thisC + d_thisC.T
            del d_thisC
            for b in range(a+1):
              Cab = astromutils.getCovar((dy,dx), PS_smooth[k,:,:], branch=(a,b))
              thisC += Cab*numpy.outer(xi[a,:],xi[b,:])
              dCdp[a,:,:] += Cab*numpy.tile(xi[b,:], (dy*dx,1))
              if b<a:
                thisC += Cab.T*numpy.outer(xi[b,:],xi[a,:])
                dCdp[b,:,:] += Cab.T*numpy.tile(xi[a,:], (dy*dx,1))
        good = numpy.where(SCAMask[k,ymin:ymax,xmin:xmax].flatten()==0)[0]
        any_good[good] = True
        good2d = numpy.ix_(good,good)
        Cinv = numpy.zeros((dy*dx,dy*dx))
        Cinv[good2d] = astromutils.syminv((CovWT[k,:,:])[good2d])

        # loop over data files
        # get datasum (to subtract)
        datasum = numpy.zeros((ny,nx))
        for ik in range(nrepeat[k]):
          data = astromutils.read_segment(re.sub(r'\*', '{:03d}'.format(ik+startInput), inputdata[k]), bounds=[ymin,ymax,xmin,xmax])\
                  /norm[k,ik]/cflat_image[ymin:ymax,xmin:xmax]
          data = data/(1.+qnorm[ty,tx]) -1.
          datasum = datasum+data
        datasum /= nrepeat[k]
        Tref = numpy.sum(Cinv*thisC) * (nrepeat[k]-1)
        myT = mysq = 0.
        #
        # now main loop
        for ik in range(nrepeat[k]):
          data = astromutils.read_segment(re.sub(r'\*', '{:03d}'.format(ik+startInput), inputdata[k]), bounds=[ymin,ymax,xmin,xmax])\
                  /norm[k,ik]/cflat_image[ymin:ymax,xmin:xmax]
          data = data/(1.+qnorm[ty,tx]) -1. - datasum
          Cid = Cinv@data.flatten()
          if ik==0: CiCthisCi = Cinv@thisC@Cinv
          for a in range(astromutils.nTemplate):
            G = Cid * (dCdp[a,:,:]@Cid)
            if ik==0: G -= (nrepeat[k]-1) * numpy.sum(dCdp[a,:,:]*CiCthisCi,axis=1)
            dNLLdp[a,:,:] += G.reshape((dy,dx))
          if ik==0: del CiCthisCi
          myT += numpy.sum(data.flatten()*Cid)
          mysq += numpy.mean(data**2)
        mysq /= nrepeat[k]-1
        #print('     T = {:12.5E} ref = {:12.5E} ratio = {:12.5E} var = {:12.5E} (th) {:12.5E} (obs)'.format(
        #  myT,Tref,myT/Tref,numpy.mean(numpy.diagonal(thisC)), mysq))

        # get response matrix (if not done already)
        cprod = numpy.zeros((astromutils.nTemplate,dy*dx,dy*dx))
        for a in range(astromutils.nTemplate):
          cprod[a,:,:] = dCdp[a,:,:]@Cinv
          for b in range(a+1):
            if a==b or not DiagFisher:
              Fisher[a*dy*dx:(a+1)*dy*dx,b*dy*dx:(b+1)*dy*dx] += (cprod[a,:,:]*cprod[b,:,:].T + Cinv*astromutils.force_sym_mult(cprod[a,:,:],dCdp[b,:,:].T))\
                * (nrepeat[k]-1) # if multiple realizations
        del cprod

        # end of k-loop

      # use symmetry for other parts of Fisher
      if not DiagFisher:
        for a in range(astromutils.nTemplate-1):
          for b in range(a+1,astromutils.nTemplate):
            Fisher[a*dy*dx:(a+1)*dy*dx,b*dy*dx:(b+1)*dy*dx] = Fisher[b*dy*dx:(b+1)*dy*dx,a*dy*dx:(a+1)*dy*dx].T

      # force flat
      if ForceFlat:
        newgrad = numpy.zeros((dy*dx,))
        any_good__ = numpy.where(any_good)
        sumall = numpy.sum((dxmap[0,ymin:ymax,xmin:xmax].flatten())[any_good__])
        newgrad[any_good__] = -sumall/flatprior**2/(dy*dx)**2
        dNLLdp[0,:,:] += newgrad.reshape((dy,dx))
        any_good2d = numpy.ix_(any_good,any_good)
        Fisher[any_good2d] += 1./flatprior**2/(dy*dx)**2

      this_change = numpy.zeros((astromutils.nTemplate, dy, dx))
      for a in range(astromutils.nTemplate):
        w = Fisher[a*dy*dx:(a+1)*dy*dx,a*dy*dx:(a+1)*dy*dx] + numpy.identity(dy*dx)/sigmaprior[a]**2
        this_change[a,:,:] = (astromutils.syminv(w)@dNLLdp[a,:,:].flatten()).reshape((dy,dx))
      #for a in range(astromutils.nTemplate):
      #  Fisher[a*dy*dx:(a+1)*dy*dx,a*dy*dx:(a+1)*dy*dx] += numpy.identity(dy*dx)/sigmaprior[a]**2
      #print(' ... start inverse', time.asctime(time.localtime(time.time())))
      #this_change = (astromutils.syminv(Fisher)@dNLLdp.flatten()).reshape((nTemplate,dy,dx))
      #print(' ...  end  inverse', time.asctime(time.localtime(time.time())))

      if cy==0: dxfiltmap[:,ymin:ymax,xmin:xmax] = this_change # report first change map
      if ChebIter:
        this_change2 = this_change/(lmax * (1. + alphaCheb*numpy.cos(numpy.pi/nCheb*((cy*kCheb)%nCheb+.5)))/2.)
      else:
        this_change2 = condition_factor*this_change
        if QuadIter:
          if cy%2==0:
            prev_change = numpy.copy(this_change)
          else:
            this_change2 += condition_gamma*prev_change

      dxmap[:,ymin:ymax,xmin:xmax] += this_change2
      print('   changes =', numpy.amax(numpy.abs(this_change2),axis=(1,2)))
      print('       rms =', numpy.sqrt(numpy.mean(this_change2**2,axis=(1,2))))
      # print temporary dxmap to a file
      hdu = fits.PrimaryHDU(dxmap.astype(numpy.float32))
      hdu.writeto(outstem+'_dx.fits', overwrite=True)
      sys.stdout.flush()
        

# initial offsets -- filtering
hdu = fits.PrimaryHDU(dxfiltmap.astype(numpy.float32))
hdu.writeto(outstem+'_dxfilt.fits', overwrite=True)
# final offsets -- filtering
hdu = fits.PrimaryHDU(dxmap.astype(numpy.float32))
hdu.writeto(outstem+'_dx.fits', overwrite=True)

