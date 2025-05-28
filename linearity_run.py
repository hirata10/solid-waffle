import sys
import re
import numpy
import time
import json
from scipy import linalg
from scipy.signal import convolve
from scipy.special import eval_legendre
from astropy.io import fits
import asdf
from datetime import datetime, timezone
from os import path

import pyirc

# get the parameters for this run
with open(sys.argv[1], 'r') as file:
    pars = json.load(file)
    script = file.read()

# get sizes of things
nside = pyirc.get_nside(pars['RAMPS'][0]['FORMAT'])
ncblock = 32 # split into blocks to speed things up -- must be a power of 2
dx = nside//ncblock
nramps = len(pars['RAMPS'])

# global sign parameter -- are we fitting for *increasing* ramps or *decreasing*?
# sign=1 for increasing
sign=-1
if 'SIGN' in pars:
   sign = int(pars['SIGN'])
print('Ramp direction:', sign)

# we'll use a lot of these
p_order = int(pars['P_ORDER'])

# is there a dark in the linearity runs?
# If so, can specify which ramp number
# (-1 for the last one)
use_dark = None
if 'DARK' in pars:
    use_dark = int(pars['DARK'])

# timiing information
tframe = 3.04
if 'TFRAME' in pars:
    tframe = float(pars['TFRAME'])

# fraction of mean slope where we say something is saturated
slopecut = .1
if 'SLOPECUT' in pars:
    slopecut = float(pars['SLOPECUT'])

# definition for getting cube
def get_median_cube(flist,fileformat=None,xmin=None,xmax=None,tstart=None,tend=None):
    # the flist is a list of files
    # tstart and tend are slices in Fortran convention
    # all the named arguments are required

    cube = numpy.zeros((tend-tstart+1,nside,dx))
    for t in range(tstart,tend+1):
        temp = numpy.zeros((nside,dx,len(flist)), dtype=numpy.float32)
        for k in range(len(flist)):
            temp[:,:,k] = pyirc.load_segment(flist[k], fileformat, [xmin,xmax,0,nside], [t], False)
        cube[t-tstart,:,:] = numpy.median(temp, axis=-1)
    if sign>0: cube = 65535-cube # flip to upward-going ramp
    return cube

# output image cube:
# 1st axis is Legendre polymonial coefficients, order 0..p_order
# then:
# p_order+1 -> Smin
# p_order+2 -> Sref
# p_order+3 -> Smax
# p_order+4 -> maxerr
# p_order+5 -> quality (0=OK)
image_all = numpy.zeros((p_order+6,nside,nside))

# pixel flats
pflat = numpy.zeros((nramps,nside,nside), dtype=numpy.float32)

# flag reference pixels
# note only the first 24 flags should be used so that this can be
# expressed exactly in IEEE 754
dq = numpy.zeros((nside,nside), dtype=numpy.uint32)
dq[:4,:] |= 2**31
dq[-4:,:] |= 2**31
dq[:,:4] |= 2**31
dq[:,-4:] |= 2**31

# get bias information
infile = pars['BIAS']['FILE']
if infile[-5:].lower()=='.asdf':
    with asdf.open(infile) as f:
        x = f
        for leaf in pars['BIAS']['PATH']:
            x = x[leaf]
        image_all[p_order+2,:,:] = x[int(pars['BIAS']['SLICE']),:,:]
else:
    raise ValueError("Only ASDF bias input currently accepted. Can't read " + infile)

# loop over strips. optional to stop somewhere
ncblock_use = ncblock
if 'STOP' in pars:
    ncblock_use=min(ncblock, int(pars['STOP']))
for i in range(ncblock_use):
    xmin = dx*i
    xmax = xmin+dx

    # get file names: flists is going to be a list of lists
    flists = []
    ntslice = numpy.zeros((nramps,), dtype=numpy.int16)
    for j in range(nramps):
        m = re.search(r'^(.+)/(\d\d\d\d\d\d\d\d)([^/]+)\_(\d+).fits', pars['RAMPS'][j]['FILE'])
        if m:
            prefix = m.group(1)
            stamp = int(m.group(2))
            body = m.group(3)
        else:
            print('Match failed.')
            exit()
        filelist = []
        for k in range(pars['RAMPS'][j]['NRAMP']):
            thisname = prefix + '/{:d}'.format(stamp) + body + '_{:03d}.fits'.format(k+pars['RAMPS'][j]['START'])
            filelist.append(thisname)
        flists.append(filelist)
        ntslice[j] = pyirc.get_num_slices(pars['RAMPS'][j]['FORMAT'], filelist[0])

    # now get the median cubes
    # also need to get the Smin and Smax
    for j in range(nramps):
        print('Strip', i, 'files -->', j, ntslice[j], flists[j]); sys.stdout.flush()
        cube = get_median_cube(flists[j],fileformat=pars['RAMPS'][j]['FORMAT'],xmin=xmin,xmax=xmax,tstart=pars['RAMPS'][j]['TSTART'],tend=ntslice[j])
        print(numpy.median(cube, axis=(1,2)))
        Smin_ = numpy.amin(cube, axis=0)
        Smax_ = numpy.amax(cube, axis=0)

        # mask later parts of the ramp if they saturate
        halfframe = (numpy.shape(cube)[0]+1)//2
        diff = (cube[halfframe,:,:]-cube[0,:,:])/halfframe
        for k in range(numpy.shape(cube)[0]-2,p_order,-1):
            #diff = (Smax_-Smin_)/(numpy.shape(cube)[0]-1) # <-- old criterion
            if sign>0:
                cuthere = cube[k+1,:,:]-cube[k-1,:,:] < 2*slopecut*diff
                Smax_ = numpy.where(cuthere, cube[k,:,:], Smax_)
            else:
                cuthere = cube[k+1,:,:]-cube[k-1,:,:] > 2*slopecut*diff
                Smin_ = numpy.where(cuthere, cube[k,:,:], Smin_)

        if j==0:
            Smin = Smin_
            Smax = Smax_
        else:
            Smin = numpy.minimum(Smin,Smin_)
            Smax = numpy.maximum(Smax,Smax_)
        print('   Smin deciles:', [numpy.round(numpy.percentile(Smin,10*ip),1) for ip in range(1,10)])
        print('   Smax deciles:', [numpy.round(numpy.percentile(Smax,10*ip),1) for ip in range(1,10)])

    # set some more quality flags
    LowSignalRange = numpy.abs(Smax-Smin)<1024
    dq[:,xmin:xmax] |= numpy.where(LowSignalRange,2**20,0).astype(numpy.uint32)

    # now extend range to cover the bias level
    if sign>0:
        Smin = numpy.minimum(Smin, image_all[p_order+2,:,xmin:xmax] - float(pars['NEGATIVEPAD']))
    else:
        Smax = numpy.maximum(Smax, image_all[p_order+2,:,xmin:xmax] + float(pars['NEGATIVEPAD']))

    # clip
    Smin = numpy.clip(Smin,0,None)
    Smax = numpy.clip(Smax,None,65535)
    # pad to avoid singularities
    Smin -= 0.5
    Smax += 0.5

    # now get the coefficients in the cost function
    # cost = sum [phi(S) - a f(S) - c]^2
    # where f(S) = timestamp
    # phi(S) = sum_L=0^p g_L P_L(z)
    # where z = S rescaled into -1 .. +1, z = (S-16383.5)/32768.
    # so for each pixel, there is a (p+1+2*nramps) x (p+1+2*nramps) matrix
    Amat = numpy.zeros((p_order+1+2*nramps,p_order+1+2*nramps,nside,dx))
    # build matrices
    for j in range(nramps):
        offset_j = p_order + 1 + 2*j # starting index for this ramp
        cube = get_median_cube(flists[j],fileformat=pars['RAMPS'][j]['FORMAT'],xmin=xmin,xmax=xmax,tstart=pars['RAMPS'][j]['TSTART'],tend=ntslice[j])
        for tt in range(numpy.shape(cube)[0]):
            z = -1. + 2.*(cube[tt,:,:]-Smin)/(Smax-Smin) # put in -1 .. +1
            mask = numpy.where(numpy.abs(z)<1,1,0) # if outside the range, mask that pixel
            u = numpy.zeros((p_order+1,nside,dx))
            for L in range(p_order+1):
                u[L,:,:] = eval_legendre(L,z)
            for L in range(p_order+1):
                for Lp in range(p_order+1):
                    Amat[L,Lp,:,:] += mask*u[L,:,:]*u[Lp,:,:]
                Amat[offset_j,L,:,:] -= mask*u[L,:,:]*tt
                Amat[L,offset_j,:,:] -= mask*u[L,:,:]*tt
                Amat[offset_j+1,L,:,:] -= mask*u[L,:,:]
                Amat[L,offset_j+1,:,:] -= mask*u[L,:,:]
            Amat[offset_j,offset_j,:,:] += mask*tt**2
            Amat[offset_j,offset_j+1,:,:] += mask*tt
            Amat[offset_j+1,offset_j,:,:] += mask*tt
            Amat[offset_j+1,offset_j+1,:,:] += mask

    # check invertibility
    D = numpy.abs(numpy.linalg.det(numpy.transpose(Amat[2:,2:,:,:],axes=(2,3,0,1))))
    Dmed = numpy.median(D)
    NearSingularMatrix = D<1e-12*Dmed
    dq[:,xmin:xmax] |= numpy.where(NearSingularMatrix,2**20,0).astype(numpy.uint32)

    # to prevent instabilities, effectively pin the higher-order coefficients in the low signal cases
    for l_ in range(2,p_order+1+2*nramps):
        Amat[l_,l_,:,:] += numpy.where(NearSingularMatrix, 1e24,0)

    # now get the coefficients assuming we start with 0, 1
    coef = numpy.zeros((p_order+1,nside,dx))
    coef[1,:,:] = 1.
    v = -numpy.linalg.solve(numpy.transpose(Amat[2:,2:,:,:],axes=(2,3,0,1)), numpy.transpose(Amat[1,2:,:,:],axes=(1,2,0)) )
    coef[2:,:,:] = numpy.transpose(v[:,:,:p_order-1], axes=(2,0,1))

    # get linearization error
    err = numpy.zeros((nside,dx))
    for j in range(nramps):
        cube = get_median_cube(flists[j],fileformat=pars['RAMPS'][j]['FORMAT'],xmin=xmin,xmax=xmax,tstart=pars['RAMPS'][j]['TSTART'],tend=ntslice[j])
        for tt in range(numpy.shape(cube)[0]):
            z = -1. + 2.*(cube[tt,:,:]-Smin)/(Smax-Smin) # put in -1 .. +1
            pred = numpy.copy(z)
            for L in range(2,p_order+1):
                pred += v[:,:,L-2]*eval_legendre(L,z)
            this_err = pred - (v[:,:,p_order-1+2*j]*tt + v[:,:,p_order+2*j])
            this_err = numpy.where(numpy.abs(z)<1, this_err, 0.)
            err = numpy.maximum(err,numpy.abs(this_err))

    # now do the rescaling to slope 1 at the bias point
    zbias = (image_all[p_order+2,:,xmin:xmax]-Smin)/(Smax-Smin)*2. - 1. # z value at the bias
    slope = numpy.zeros((nside,dx))
    for L in range(1,p_order+1):
        # this is to get dP_L(z)/dz, which is equal to sum_{L'<L, L-L' odd} (2L'+1) P_L'(z)
        dPl = numpy.zeros((nside,dx))
        for Lp in range((L-1)%2,L,2):
            dPl += (2*Lp+1)*eval_legendre(Lp,zbias)
        slope += coef[L,:,:] * dPl
        del dPl
    slope *= 2./(Smax-Smin) # this is dz/dS
    coef = coef/slope[None,:,:]
    err = err/slope # also need to re-scale the fit error - this is now in linearized DN

    print('   err deciles:', [numpy.round(numpy.percentile(err,10*ip),1) for ip in range(1,10)])

    # get the P-flat fields
    for j in range(nramps):
        pflat[j,:,xmin:xmax] = (v[:,:,p_order-1+2*j]/slope/tframe).astype(numpy.float32)

    del slope

    # and now the intercept
    icpt = numpy.zeros((nside,dx))
    for L in range(1,p_order+1):
        icpt += coef[L,:,:] * eval_legendre(L,zbias)
    coef[0,:,:] = -icpt
    del icpt

    # monotonicity checks
    nzstep = 512 # number of grid points to check monotonicity
    for zstep in range(nzstep+1):
        z = -1. + 2*zstep/nzstep
        Sprime = numpy.zeros((nside,dx))
        for L in range(1,p_order+1):
            # this is to get dP_L(z)/dz, which is equal to sum_{L'<L, L-L' odd} (2L'+1) P_L'(z)
            dPl = numpy.zeros((nside,dx))
            for Lp in range((L-1)%2,L,2):
                dPl += (2*Lp+1)*eval_legendre(Lp,zbias)
            Sprime += coef[L,:,:] * dPl
            del dPl
        dq[:,xmin:xmax] |= numpy.where(Sprime<=0,2**20,0).astype(numpy.uint32)

    image_all[:p_order+1,:,xmin:xmax] = coef
    image_all[p_order+1,:,xmin:xmax] = Smin
    image_all[p_order+3,:,xmin:xmax] = Smax
    image_all[p_order+4,:,xmin:xmax] = err

# 4D or 3D to 2D projection, displaying lots of strips left to right for visualization
# (not actually needed in the final script)
def to2D(arr):
    d = numpy.shape(arr)
    if len(d)==4:
        return numpy.transpose(arr, axes=(2,0,1,3)).reshape(d[2],-1)
    else:
        return numpy.transpose(arr, axes=(1,0,2)).reshape(d[1],-1)

# now OK to convert to float32
image_all = image_all.astype(numpy.float32)

# quality flags store exactly in IEEE754 32-bit arithmetic
image_all[p_order+5,:,:] = dq.astype(numpy.float32)

# save the data cube
Im = fits.PrimaryHDU(image_all)
Im.header['P_ORDER'] = p_order
Im.header['NSIDE_Y'] = nside
Im.header['NSIDE_X'] = nside
for L in range(p_order):
    Im.header['SLICE{:03d}'.format(L+1)] = 'Legendre coef L={:d}'.format(L)
Im.header['SLICE{:03d}'.format(p_order+1)] = 'min signal in fit'
Im.header['SLICE{:03d}'.format(p_order+2)] = 'ref signal in fit (usually bias level; slope=1)'
Im.header['SLICE{:03d}'.format(p_order+3)] = 'max signal in fit'
Im.header['SLICE{:03d}'.format(p_order+4)] = 'max error'
Im.header['SLICE{:03d}'.format(p_order+5)] = 'quality'
Im.header['QBIT20'] = 'Bit 20: failed or unstable nonlinearity solution'
Im.header['QBIT31'] = 'Bit 31: ref pixel'

# do the dark subtraction
dark = numpy.zeros((nside,nside), dtype=numpy.float32)
if use_dark is not None:
    dark = pflat[use_dark,:,:]
    pflat = pflat - dark[None,:,:]

# quality flag information
for j in range(32):
    print('flag {:2d}, count {:7d}'.format(j, numpy.count_nonzero(dq&(2**j))))

fits.HDUList([Im, fits.ImageHDU(pflat)]).writeto(pars['OUTPUT'][:-5]+'.fits', overwrite=True)

# where the input data came from ...
pedigree = 'DUMMY'
if 'PEDIGREE' in pars:
    pedigree = pars['PEDIGREE']

# ASDF output
tree = {'roman': {
    'meta': {
        'author': 'linearity_run.py',
        'description': 'linearity_run.py',
        'instrument': {
            'detector': 'WFI{:02d}'.format(int(pars['SCA'])),
            'name': 'WFI'
        },
        'origin': 'PIT - solid-waffle',
        'date': datetime.now(timezone.utc).isoformat(),
        'pedigree': pedigree,
        'reftype': 'LINEARITYLEGENDRE',
        'telescope': 'ROMAN',
        'useafter': '!time/time-1.2.0 2020-01-01T00:00:00.000'
    },
    'data': image_all[:p_order+1,:,:],
    'pflat': pflat,
    'dark': dark,
    'dq': dq,
    'Smin': image_all[p_order+1,:,:],
    'Sref': image_all[p_order+2,:,:],
    'Smax': image_all[p_order+3,:,:],
    'maxerr': image_all[p_order+4,:,:]
},
'notes': {
    'script': script,
    'units': "data, maxerr: DN_lin; pflat, dark: DN_lin/s; Smin, Sref, Smax: DN_raw"
}
}

asdf.AsdfFile(tree).write_to(pars['OUTPUT'][:-5]+'.asdf')
