import sys
import numpy
import scipy
import astropy
from astropy.io import fits
import scipy.stats
import scipy.ndimage
import fitsio
import copy
import warnings
from fitsio import FITS,FITSHDR
from ftsolve import center,decenter,solve_corr,solve_corr_many,solve_corr_vis,solve_corr_vis_many,pad_to_N
from scipy.signal import correlate2d,fftconvolve

# <== TESTING PARAMETERS ONLY ==>
#
# [these are false and should only be set to true for debugging purposes]
Test_SubBeta = False

# <== THESE FUNCTIONS DEPEND ON THE FORMAT OF THE INPUT FILES ==>

# Version number of script
def get_version():
  return 27

# Function to get array size from format codes in load_segment
# (Note: for WFIRST this will be 4096, but we want the capability to
# run this script on H1/H2RG data.)
#
def get_nside(formatpars):
  if formatpars==1: return 4096
  if formatpars==2: return 2048
  if formatpars==3: return 4096
  if formatpars==4: return 4096

# Get number of time slices
def get_num_slices(formatpars, filename):

  # Switch based on input format
  if formatpars==1 or formatpars==2:
    hdus = fits.open(filename)
    ntslice = int(hdus[0].header['NAXIS3'])
    hdus.close()
  elif formatpars==3:
    hdus = fits.open(filename)
    ntslice = len(hdus)-1
    hdus.close()
  elif formatpars==4:
    hdus = fits.open(filename)
    ntslice = int(hdus[1].header['NAXIS3'])
    hdus.close()
  else:
    print ('Error! Invalid formatpars =', formatpars)
    exit()
  return ntslice

# Function to load an image segment
#
# filename = name of the source FITS file
# formatpars = integer describing which type of format to use
#     format 1: H4RG, all data in a single HDU, ramp slope positive (ex. DCL H4RG-18237 data)
# xyrange = list [xmin,xmax,ymin,ymax] (first row/col are zero) -- EXcluding xmax and ymax!
# tslices = list of time slices to use (first is *1*)
# verbose = True or False (use True only for de-bugging)
#
# Returns a 3D array of dimension number tslices, ymax-ymin, xmax-xmin
#
def load_segment(filename, formatpars, xyrange, tslices, verbose):
  if verbose: print ('Reading:', filename)

  # Recommended True (False defaults to astropy tools, which work but are slow because of the way this script works)
  use_fitsio = True

  # Get dimensions of output cube
  nxuse = xyrange[1]-xyrange[0]
  nyuse = xyrange[3]-xyrange[2]
  ntslice_use = len(tslices)
  output_cube = numpy.zeros((ntslice_use, nyuse, nxuse))

  # Switch based on input format
  if formatpars==1 or formatpars==2:
    if use_fitsio:
      fileh = fitsio.FITS(filename)
      N = get_nside(formatpars)
      for ts in range(ntslice_use):
        t = tslices[ts]
        if ts>0 and tslices[ts]==tslices[ts-1]:
          output_cube[ts,:,:] = output_cube[ts-1,:,:] # don't read this slice again
        else:
          output_cube[ts,:,:] = 65535 - numpy.array(fileh[0][t-1, xyrange[2]:xyrange[3], xyrange[0]:xyrange[1]])
      fileh.close()
    else:
      hdus = fits.open(filename)
      in_hdu = hdus[0]
      ntslice = in_hdu.data.shape[0]
      if verbose:
        print ('input shape -> ', in_hdu.data.shape)
        print ('number of slices =', ntslice, ', used =', ntslice_use)
      for ts in range(ntslice_use):
        t = tslices[ts]
        output_cube[ts,:,:] = 65535 - in_hdu.data[t-1, xyrange[2]:xyrange[3], xyrange[0]:xyrange[1]]
      hdus.close()
  elif formatpars==3:
    if use_fitsio:
      fileh = fitsio.FITS(filename)
      N = get_nside(formatpars)
      for ts in range(ntslice_use):
        t = tslices[ts]
        output_cube[ts,:,:] = numpy.array(fileh[t][xyrange[2]:xyrange[3], xyrange[0]:xyrange[1]])
      fileh.close()
    else:
      print ('Error: non-fitsio methods not yet supported for formatpars=3')
      exit()
  elif formatpars==4:
    if use_fitsio:
      fileh = fitsio.FITS(filename)
      N = get_nside(formatpars)
      for ts in range(ntslice_use):
        t = tslices[ts]
        output_cube[ts,:,:] = numpy.array(fileh[1][0, t-1, xyrange[2]:xyrange[3], xyrange[0]:xyrange[1]])
      fileh.close()
    else:
      print ('Error: non-fitsio methods not yet supported for formatpars=4')
      exit()
  else:
    print ('Error! Invalid formatpars =', formatpars)
    exit()

  return output_cube

# <== FUNCTIONS BELOW HERE ARE INDEPENDENT OF THE INPUT FORMAT ==>

# Dictionary of indices
# These are designed for consistency with the outputs lists of certain functions
class IndexDictionary:
  def __init__(self, itype):
    # basic characterization parameter index list -- outputs for "basic" function
    if itype==0:
      self.Nb = 11 # number of basic parameters
      #
      self.ngood = 0
      self.graw = 1
      self.gacorr = 2
      self.g = 3
      self.alphaH = 4
      self.alphaV = 5
      self.beta = 6
      self.I = 7
      self.alphaD = 8
      self.tCH = 9
      self.tCV = 10
      #
      # polychar indices: output is list [1, ***, residual]
      # intermediate terms are basic output [ind1:ind2]
      self.ind1 = 3
      self.ind2 = 9
      self.indp1 = 1
      self.indp2 = self.indp1 + self.ind2-self.ind1
      #
      self.N = self.Nbb = self.Nb

  # adds BFE kernel, (2s+1) x (2s+1)
  def addbfe(self, s):
    self.s = s
    self.ker0 = self.Nb + 2*s*(s+1)  # BFE (0,0) kernel index
    self.Nbb += (2*s+1)*(2*s+1)
    self.N += (2*s+1)*(2*s+1)

  # adds higher-order non-linearity coefficients (p coefs, for total degree 1+p)
  def addhnl(self, p):
    self.p = p
    self.N += p

swi = IndexDictionary(0)

# Routine to get percentile cuts with a mask removed
#
# mask consists of 0's and 1's and is the same size as this_array
def pyIRC_percentile(this_array, mask, perc):
  val = this_array.flatten()
  ma = mask.flatten()
  w = numpy.array([val[x] for x in numpy.where(ma>.5)])
  return numpy.percentile(w,perc)

# Routine to get mean with a mask removed
def pyIRC_mean(this_array, mask):
  val = this_array.flatten()
  ma = mask.flatten()
  w = numpy.array([val[x] for x in numpy.where(ma>.5)])
  return numpy.mean(w)

# Get reference corrections from left & right pixel sets
# yrange = [ymin,ymax] (inclusive)
#
# Output depends on the length of tslices:
#  elements 0 .. ntslice_use-1 -> median of that time slice
#  elements ntslice_use .. 2*ntslice_use-1 -> median of (first) - (this slice)
#  (if ntslice_use>=2) then
#    element 2*ntslice_use -> median of [(-2) - (-1)] - [(0) - (1)] (otherwise this is 0)
#    [this is really used to measure curvature of the reference pixel ramp]
#
# output always has length 2*ntslice_use+1
#
def ref_corr(filename, formatpars, yrange, tslices, verbose):

  # Side length of the array (needed to find reference pixel indexing)
  N = get_nside(formatpars)
  # Number of time slices
  ntslice_use = len(tslices)
  # Clear list
  output_ref = []

  # Build arrays of reference pixels
  my_array_L = load_segment(filename, formatpars, [0,4]+yrange, tslices, False)
  my_array_R = load_segment(filename, formatpars, [N-4,N]+yrange, tslices, False)
  my_array_LR = numpy.concatenate((my_array_L, my_array_R), axis=2)
  if verbose: print (N, my_array_LR.shape)

  for ts in range(ntslice_use):
    output_ref.append(numpy.median(my_array_LR[ts,:,:]))
  for ts in range(ntslice_use):
    diff_array = my_array_LR[0,:,:] - my_array_LR[ts,:,:]
    output_ref.append(numpy.median(diff_array))
  if ntslice_use>1:
    diff_array = my_array_LR[ntslice_use-2,:,:] - my_array_LR[ntslice_use-1,:,:]\
                 -(my_array_LR[0,:,:]-my_array_LR[1,:,:])*(tslices[-1]-tslices[-2])/float(tslices[1]-tslices[0])
    output_ref.append(numpy.median(diff_array))
  else:
    output_ref.append(0)
  return output_ref

#
# Get reference corrections from left & right pixel sets
# for a full list of files.
# ny = number of y-bins (e.g. 32 for an H4RG and regions of 128 pixel size in y-direction)
#
# Output depends on the length of tslices:
#  elements 0 .. ntslice_use-1 -> median of that time slice
#  elements ntslice_use .. 2*ntslice_use-1 -> median of (first) - (this slice)
#  (if ntslice_use>=2) then
#    element 2*ntslice_use -> median of [(-2) - (-1)] - [(0) - (1)] (otherwise this is 0)
#    [this is really used to measure curvature of the reference pixel ramp]
#
# output is stored in a numpy array of size num_files, ny, 2*ntslice_use+1
#
def ref_array(filelist, formatpars, ny, tslices, verbose):

  num_files = len(filelist)
  ntslice_use = len(tslices)
  output_array = numpy.zeros((num_files, ny, 2*ntslice_use+1))

  dy = get_nside(formatpars)//ny
  for ifile in range(num_files):
    for iy in range(ny):
      ymin = dy*iy
      ymax = ymin+dy
      output_array[ifile, iy, :] = numpy.asarray(ref_corr(filelist[ifile], formatpars, [ymin,ymax], tslices, False))
      if verbose:
        print (ifile, iy)
        print (output_array[ifile, iy, :])

  return(output_array)
#
# Similar but if we only need one row (iy) to be good
# *** Only use this function if you are absolutely sure of what you need!
def ref_array_onerow(filelist, formatpars, iy, ny, tslices, verbose):
  num_files = len(filelist)
  ntslice_use = len(tslices)
  output_array = numpy.zeros((num_files, ny, 2*ntslice_use+1))
  dy = get_nside(formatpars)//ny
  for ifile in range(num_files):
    ymin = dy*iy
    ymax = ymin+dy
    output_array[ifile, iy, :] = numpy.asarray(ref_corr(filelist[ifile], formatpars, [ymin,ymax], tslices, False))
    if verbose:
      print (ifile, iy)
      print (output_array[ifile, iy, :])
  return(output_array)
#
# similar but uses a user-specified range of y-values, and output lacks the 'iy' index
# (i.e. is 2D array)
def ref_array_block(filelist, formatpars, yrange, tslices, verbose):

  num_files = len(filelist)
  ntslice_use = len(tslices)
  output_array = numpy.zeros((num_files, 2*ntslice_use+1))

  if len(yrange)<2:
    print ('Error in ref_array_block: yrange =', yrange)
    exit()
  for ifile in range(num_files):
    ymin = yrange[0]
    ymax = yrange[1]
    output_array[ifile, :] = numpy.asarray(ref_corr(filelist[ifile], formatpars, [ymin,ymax], tslices, False))
    if verbose:
      print (ifile)
      print (output_array[ifile, :])

  return(output_array)

# Generate a 4D date cube containing information on a region of the detector
#
# filename = name of the source FITS file
# formatpars = integer describing which type of format to use
#     format 1: H4RG, all data in a single HDU, ramp slope positive (ex. DCL H4RG-18237 data)
# xyrange = list [xmin,xmax,ymin,ymax] (first row/col are zero) -- EXcluding xmax and ymax!
# tslices = list of time slices to use (first is *1*)
# maskinfo = information on how the masking works (list format, if not enough elements goes to default)
#   maskinfo[0] = range around median to accept (default: 0.1, must be within 10% of median)
#   maskinfo[1] = boolean, mask assuming light exposure (default: True)
#
# verbose = True or False (use True only for de-bugging)
#
# Returns a 4D array of dimension number of files +1, number tslices, ymax-ymin, xmax-xmin
#   the *last* "file" is the mask (0 or 1)
#
def pixel_data(filelist, formatpars, xyrange, tslices, maskinfo, verbose):

  # Masking parameters
  cut_offset = 0.1
  if len(maskinfo)>=1: cut_offset = maskinfo[0]
  do_mask = True
  if len(maskinfo)>=2: do_mask = maskinfo[1]

  num_files = len(filelist)
  ntslice_use = len(tslices)
  output_array = numpy.zeros((num_files+1, ntslice_use, xyrange[3]-xyrange[2], xyrange[1]-xyrange[0]))

  for ifile in range(num_files):
    output_array[ifile,:,:,:] = load_segment(filelist[ifile], formatpars, xyrange, tslices, verbose)

  # Generate mean CDS image and consider the median
  mCDS = numpy.mean(output_array[0:num_files,0,:,:], axis=0) - numpy.mean(output_array[0:num_files,-1,:,:], axis=0)
  mCDS_med = numpy.median(mCDS)
  if do_mask:
    a = (1./mCDS_med)*mCDS
    goodmap = numpy.where(numpy.logical_and(a>1-cut_offset,a<1+cut_offset),1,0)
  else:
    goodmap = numpy.ones_like(mCDS)
  for f in range(num_files):
    for t in range(ntslice_use):
      goodmap *= numpy.where(output_array[f,t,:,:]>0,1,0)
  if verbose:
    print ('Median =', mCDS_med, 'cut_offset =', cut_offset)
    print (goodmap)
    print (goodmap.shape)
  # Copy map of good pixels into the output
  for t in range(ntslice_use):
    output_array[num_files,t,:,:] = goodmap

  return output_array

# Routine to get nonlinearity curve
# Inputs:
#   filelist       <- list of 'light' files
#   formatpars     <- format type for file
#   timeslice      <- integer (does frames 1 .. tmax) or list [tref,t1,t2], uses t1 ... t2 with reset at tref
#   ngrid          <- number of cells, list [ny,nx]
#   Ib             <- shouldn't need this except for de-bugging (forces a specific quadratic fit)
#   usemode        <- 'dev' (deviation from beta fit) or 'abs' (absolute -- zero of time is absolute)
#   verbose        <- Boolean
#
# Output:
#   signal, fit, derivative [, coefs]
#   array of dimensions [tmax, ny, nx] containing reference corrected signal (in DN)
#     This is median within a file, followed by the mean.
#   The second and third outputs have the same shape as the first:
#     2nd is a polynomial fit
#     3rd is a corrected derivative (in DN/frame)
#   The polynomial coefficients are reported as coefs (ascending order: t^0, then t^1 ...) in 'abs' mode,
#     and as a data cube [order+1, ny, nx]
def gen_nl_cube(filelist, formatpars, timeslice, ngrid, Ib, usemode, verbose):
  # Extract basic information
  nfiles = len(filelist)
  nx = ngrid[1]; ny = ngrid[0]
  N = get_nside(formatpars)
  dx = N//nx; dy = N//ny

  # Check whether we have a list or single slice
  if isinstance(timeslice, list):
    tref = timeslice[0]
    tmin = timeslice[1]
    tmax = timeslice[2]
    nt = tmax-tmin+1
  else:
    tref = 0
    tmin = 1
    nt = tmax = timeslice

  output_array = numpy.zeros((nt, ny, nx))
  temp_array = numpy.zeros((nfiles, ny, nx))

  # order of polynomial fit per pixel
  my_order = 5
  if usemode=='abs': my_order = swi.p

  if verbose:
    print ('Nonlinear cube:')
    sys.stdout.write('  reference pixel extraction ...'); sys.stdout.flush()

  # Extract reference information
  # Now ref_signal[ifile, iy, it] contains it_th time slice of group iy of ref pixels in file ifile
  ref_signal = ref_array(filelist, formatpars, ny, range(tmin,tmax+1), False)

  if verbose:
    print ('  done.')
    sys.stdout.write('Time slices:'); sys.stdout.flush()

  # Now loop over times
  for t in range(tmin,tmax+1):
    temp_array[:,:,:] = 0.
    if verbose:
      sys.stdout.write(' {:2d}'.format(t)); sys.stdout.flush()
    for ifile in range(nfiles):
      val = load_segment(filelist[ifile], formatpars, [0,N,0,N], [1,t], False) # make 2D array
      valc = val[1,:,:] - val[0,:,:]
      for iy in range(ny):
        for ix in range(nx):
          temp_array[ifile,iy,ix] = numpy.median(valc[dy*iy:dy*(iy+1), dx*ix:dx*(ix+1)])\
            - (ref_signal[ifile, iy, t-tmin] - ref_signal[ifile, iy, 0])
    output_array[t-tmin,:,:] = -numpy.mean(temp_array,axis=0)
    # <-- note: we flipped the sign so that the signal is positive

  if verbose: print ('')

  # Make fit and derivatives
  coefs_array = numpy.zeros((my_order+1, ny, nx))
  fit_array = numpy.zeros_like(output_array)
  deriv_array = numpy.zeros_like(output_array)
  for iy in range(ny):
    for ix in range(nx):
      p = numpy.poly1d(numpy.polyfit(numpy.asarray(range(tmin-tref,tmax+1-tref)), output_array[:,iy,ix], my_order))
      # p = numpy.poly1d([-Ib[iy,ix],1,0]) # <-- force proportional to beta-model (if not commented; for debugging only)
      q=numpy.poly1d.deriv(p)
      fit_array[:,iy,ix] = p(range(tmin-tref,tmax+1-tref))
      deriv_array[:,iy,ix] = q(range(tmin-tref,tmax+1-tref))
      coefs_array[:p.order+1,iy,ix] = p.c[::-1]
  if usemode=='dev':
    return output_array, fit_array, deriv_array
  else:
    return output_array, fit_array, deriv_array, coefs_array

# Routine to estimate the fractional gain error,
# log( gain[full NL] / gain[est. quad.])
# caused by using a beta-model for the nonlinearity curve instead of the full curve.
#
# Inputs are:
#   fit_array = signal in DN for true curve (length tmax array, starting with frame #1)
#   deriv_array = signal rate in DN/frame
#   Ib = charge per frame times beta (unitless)
#   tslices = time slices used for the gain (select 3 here)
#   reset_frame = which frame was reset? (0 if 1st frame is 1 after reset)
#
def compute_gain_corr(fit_array, deriv_array, Ib, tslices, reset_frame):
  # unpack time information
  ta = tslices[0] - reset_frame
  tb = tslices[1] - reset_frame
  td = tslices[2] - reset_frame
  # indices
  ina = tslices[0]-1
  inb = tslices[1]-1
  ind = tslices[2]-1

  # We want the nonlinearity corrections (mean[td]-mean[tb])/(var[tad]-var[tab])
  # which, for a non-linearity curve f(t) with f(0)=0, f'(0)=1, is
  # [f(td)-f(tb)] / { [ta*(f'(td)-f'(ta))^2 + tad*f'(td)^2] - [ta*(f'(tb)-f'(ta))^2 + tab*f'(tb)^2] }
  # = [f(td)-f(tb)] / [ ta*(f'(td)^2-f'(tb)^2 - 2f'(ta)*(f'(td)-f'(tb)) ) + td*f'(td)^2 - ta*f'(td)^2 - tb*f'(tb)^2 + ta*f'(tb)^2 ]
  # = [f(td)-f(tb)] / [ -2 ta f'(ta) (f'(td)-f'(tb)) + td f'(td)^2 - tb f'(tb)^2 ]
  #
  # Let us call that factor e^{epsilon}
  # Now for the beta-model, f(t) = t - Ib * t^2
  # so e^{epsilon} = tbd [1 - Ib (tb+td) ] / [ 4 Ib ta tbd (1 - 2 Ib ta) + td (1 - 2 Ib td)^2 - tb (1 - 2 Ib tb)^2 ]
  #                = [1 - Ib (tb+td) ] / [ 4 Ib ta (1 - 2 Ib ta) + 1 - 4 Ib (tb+td) + Ib^2 (td^2 + tb td + tb^2) ]
  # to lowest order in Ib (what is in the notes):
  # epsilon ~ Ib (-4ta+3tb+3td)

  true_expepsilon = (fit_array[ind]-fit_array[inb]) / (-2*ta*deriv_array[ina]*(deriv_array[ind]-deriv_array[inb])
    + td*deriv_array[ind]**2 - tb*deriv_array[inb]**2)
  return numpy.log(true_expepsilon*deriv_array[0]) - Ib*(-4.*ta+3.*tb+3.*td)
#
# Same but for a whole grid of points (i.e. all nx x ny cells at a time).
# also includes a mask to avoid error messages.
def compute_gain_corr_many(fit_array, deriv_array, Ib, tslices, reset_frame, is_good):
  out_array = numpy.zeros_like(fit_array[0,:,:])
  ny = numpy.shape(fit_array)[1]; nx = numpy.shape(fit_array)[2]
  for iy in range(ny):
    for ix in range(nx):
      if is_good[iy,ix]>.5:
        out_array[iy,ix] = compute_gain_corr(fit_array[:,iy,ix], deriv_array[:,iy,ix], Ib[iy,ix], tslices, reset_frame)
  return out_array

# Routine to estimate the correction to the correlation
#   g^2 C_{abab}(+/-1,0)/(It_{ab}) / (2 alpha (1-4 alpha) )
# caused by using a beta-model for the nonlinearity curve instead of the full curve.
#
# Inputs are:
#   fit_array = signal in DN for true curve (length tmax array, starting with frame #1)
#   deriv_array = signal rate in DN/frame
#   Ib = charge per frame times beta (unitless)
#   tslices = time slices used for the x-corr (select 2 here)
#   reset_frame = which frame was reset? (0 if 1st frame is 1 after reset)
#
def compute_xc_corr(fit_array, deriv_array, Ib, tslices, reset_frame):
  # unpack time information
  ta = tslices[0] - reset_frame
  tb = tslices[1] - reset_frame
  # indices
  ina = tslices[0]-1
  inb = tslices[1]-1

  # We want the correction
  # f'(tb)^2 + ta/tab * (f'(tb)-f'(ta)) - (1 - 4 Ib tb)
  return( (deriv_array[inb]**2 - ta/(tb-ta)*(deriv_array[inb]-deriv_array[ina])**2) / deriv_array[0]**2
    - (1. - 4*Ib*tb) )
#
# Same but for a whole grid of points (i.e. all nx x ny cells at a time).
# also includes a mask to avoid error messages.
def compute_xc_corr_many(fit_array, deriv_array, Ib, tslices, reset_frame, is_good):
  out_array = numpy.zeros_like(fit_array[0,:,:])
  ny = numpy.shape(fit_array)[1]; nx = numpy.shape(fit_array)[2]
  for iy in range(ny):
    for ix in range(nx):
      if is_good[iy,ix]>.5:
        out_array[iy,ix] = compute_xc_corr(fit_array[:,iy,ix], deriv_array[:,iy,ix], Ib[iy,ix], tslices, reset_frame)
  return out_array

# Routine to get IPC-corrected gain
# 
# Inputs:
#   graw      = uncorrected gain (e/DN)
#   CH        = horizontal correlation (DN^2)
#   CV        = vertical correlation (DN^2)
#   signal    = signal in this ramp (DN)
#
# Output list:
#   gain (alpha corr), e/DN
#   alphaH
#   alphaV
#
# returns [] if failed.
def gain_alphacorr(graw, CH, CV, signal):
  g = graw
  for i in range(100):
    alphaH = CH*g/(2*signal)
    alphaV = CV*g/(2*signal)
    if (alphaH+alphaV>0.25): return [] # FAIL!
    g = graw*( (1-2*(alphaH+alphaV))**2 + 2*(alphaH**2+alphaV**2) )
  return [g, alphaH, alphaV]

# Routine to get IPC+NL-corrected gain
#
# Inputs:
#   graw        = uncorrected gain (e/DN)
#   CH          = horizontal correlation (DN^2)
#   CV          = vertical correlation (DN^2)
#   signal      = signal in this ramp (DN)
#   frac_dslope = mean signal rate in (cd) / mean signal rate in (ab) - 1
#   times       = list of times [a,b,c,d] used, normalized to reference slice
#
# Output list:
#   gain g (alpha corr), e/DN
#   alphaH
#   alphaV
#   beta
#   current I (electrons per time slice)
#
# returns [] if failed
def gain_alphabetacorr(graw, CH, CV, signal, frac_dslope, times):

  # This is solving the following set of equations
  # (see Hirata's brighter-fatter effect paper)
  #
  # graw = g * [ 1 + beta I (3tb+3td-4ta) ] / [ (1-4alpha)^2 + 2alphaH^2 + 2alphaV^2 ]
  # CH = (2 I tad alphaH / g^2) [ 1 - 4alpha - 4 beta I td ]
  # CV = (2 I tad alphaV / g^2) [ 1 - 4alpha - 4 beta I td ]
  # signal = I tad [ 1 - beta I (ta+td) ] / g
  # frac_dslope = - beta I (tc+td-ta-tb)

  # Initial guess
  g = graw
  alpha = alphaH = alphaV = beta = 0
  I = signal*g/(times[3]-times[0])

  # Iterate
  # (100 iterations is overkill for this problem if alpha and beta are small)
  for numIter in range(100):
    g = graw * ((1-4*alpha)**2+2*(alphaH**2+alphaV**2)) / (1+beta*I*(3*(times[1]+times[3])-4*times[0]))
    if g<1e-3:
      print ('Gain did not converge')
      print ('IN:', graw, CH, CV, signal, frac_dslope, times)
      print ('STATUS:', g, alphaH, alphaV, alpha, I, beta)
      exit()
    temp = (1-4*alpha-4*beta*I*times[3])*2*I*(times[3]-times[0])/g**2
    alphaH = CH/temp
    alphaV = CV/temp
    if (alphaH+alphaV>0.25): return [] # FAIL!
    alpha = (alphaH+alphaV)/2.
    I = signal*g/(times[3]-times[0])/(1-beta*I*(times[3]+times[0]))
    beta = -frac_dslope/I/(times[2]+times[3]-times[0]-times[1])
    if numpy.fabs(beta)*I*(times[3]+times[0])>0.5: return [] # FAIL!

  return [g, alphaH, alphaV, beta, I]

# Basic characterization of a data cube
#
# region_cube = 4D array of the region of interest (order of indices: file, timeslice, y-ymin, x-xmin)
# dark_cube = same, but for a suite of darks
# tslices = list of time slices
# lightref = reference pixel table for correcting light exposures (2D)
#   size is num_files, 2*ntslice_use+1 (assumes we are taking the correct y-slice)
# darkref = same as for lightref
# ctrl_pars = control parameter list
#   ctrl_pars.epsilon = cut fraction (default to 0.01)
#   ctrl_pars.subtr_corr = mean subtraction for the IPC correlation? (default to True)
#   ctrl_pars.noise_corr = noise subtraction for the IPC correlation? (default to True)
#   ctrl_pars.reset_frame = reset frame (default to 0)
#   ctrl_pars.subtr_href = reference pixel subtraction? (default to True)
#   ctrl_pars.full_corr = which parameters to report (default to True = standard basic pars; False = correlation data instead)
#   ctrl_pars.leadtrailSub = lead-trail subtraction? (default to False)
#   ctrl_pars.g_ptile = percentile for inter-quantile range (default to 75)
# verbose = True or False  (recommend True only for de-bugging)
#
# Returns a list of basic calibration parameters.
# if ctrl_pars.full_corr is True
#   [number of good pixels, gain_raw, gain_acorr, gain_abcorr, aH, aV, beta, I, 0., tCH, tCV]
# if False:
#   [number of good pixels, median, variance, tCH, tCV, tCD]
# Returns the null list [] if failed.
#
# Includes a test so this won't crash if tslices[1]>=tslices[-1] but returns meaningful x-correlation C_{abab}
# (everything else is nonsense in this case)
#
def basic(region_cube, dark_cube, tslices, lightref, darkref, ctrl_pars, verbose):

  # Settings:
  newMeanSubMethod = True     # use False only for test/debug
  leadtrailSub = True         # subtract leading & trailing (by +/-4 pix) from horiz & vert correlations

  g_ptile = 75.               # percentile use for inter-quantile range for variance (default: 75, giving standard IQR)

  # Extract basic parameters
  num_files = region_cube.shape[0]-1
  nt = region_cube.shape[1]
  dy = region_cube.shape[2]
  dx = region_cube.shape[3]
  npix = dx*dy
  if nt!=len(tslices):
    print ('Error in pyirc.basic: incompatible number of time slices')
    exit()
  if verbose: print ('nfiles = ',num_files,', ntimes = ',nt,', dx,dy=',dx,dy)
  treset = 0
  if hasattr(ctrl_pars,'reset_frame'): treset = ctrl_pars.reset_frame

  # First get correlation parameters
  epsilon = .01
  if hasattr(ctrl_pars,'epsilon'): epsilon = ctrl_pars.epsilon
  subtr_corr = True
  if hasattr(ctrl_pars,'subtr_corr'): subtr_corr = ctrl_pars.subtr_corr
  noise_corr = True
  if hasattr(ctrl_pars,'noise_corr'): noise_corr = ctrl_pars.noise_corr
  if verbose: print ('corr pars =', epsilon, subtr_corr, noise_corr)
  #

  # Reference pixel subtraction?
  subtr_href = True
  if hasattr(ctrl_pars,'subtr_href'): subtr_href = ctrl_pars.subtr_href

  # return full correlation information?
  full_corr = True
  if hasattr(ctrl_pars,'full_corr'): full_corr = ctrl_pars.full_corr

  # lead-trail subtraction for IPC correlations?
  if hasattr(ctrl_pars,'leadtrailSub'): leadtrailSub = ctrl_pars.leadtrailSub

  # quantile for variance?
  if hasattr(ctrl_pars,'g_ptile'): g_ptile = ctrl_pars.g_ptile

  # Get means and variances at the early and last slices
  # (i.e. 1-point information)
  gauss_iqr_in_sigmas = scipy.stats.norm.ppf(g_ptile/100.)*2  # about 1.349 for g_ptile=75.
  box1 = region_cube[0:num_files,0,:,:] - region_cube[0:num_files,1,:,:]
  box2 = region_cube[0:num_files,0,:,:] - region_cube[0:num_files,-1,:,:]
  box2Noise = dark_cube[0:num_files,0,:,:] - dark_cube[0:num_files,-1,:,:]
  #
  if subtr_href:
    for f in range(num_files):
      if verbose: print ('lightref.shape=',lightref.shape, 'subtr ->', lightref[f,nt+1], lightref[f,2*nt-1], darkref[f,2*nt-1])
      box1[f,:,:] -= lightref[f,nt+1]
      box2[f,:,:] -= lightref[f,2*nt-1]
      box2Noise[f,:,:] -= darkref[f,2*nt-1]
  mean1 = numpy.mean(box1, axis=0)
  mean2 = numpy.mean(box2, axis=0)
  med1 = numpy.median(mean1)
  med2 = numpy.median(mean2)
  var1 = 0
  var2 = 0
  corr_mask = region_cube[-1,0,:,:]
  for if1 in range(1,num_files):
    for if2 in range(if1):
      temp_box = box1[if1,:,:] - box1[if2,:,:]
      iqr1 = pyIRC_percentile(temp_box,corr_mask,g_ptile) - pyIRC_percentile(temp_box,corr_mask,100-g_ptile)
      temp_box = box2[if1,:,:] - box2[if2,:,:]
      iqr2 = pyIRC_percentile(temp_box,corr_mask,g_ptile) - pyIRC_percentile(temp_box,corr_mask,100-g_ptile)
      var1 += (iqr1/gauss_iqr_in_sigmas)**2/2.
      var2 += (iqr2/gauss_iqr_in_sigmas)**2/2.
      if verbose: print ('Inner loop,', if1, if2, temp_box.shape)
  var1 /= num_files*(num_files-1)/2.
  var2 /= num_files*(num_files-1)/2.
  if var2<=var1 and tslices[1]<tslices[-1]: return [] # FAIL!
  gain_raw = (med2-med1)/(var2-var1+1e-100) # in e/DN
    # 1e-100 does nothing except to prevent an error when var1 and var2 are exactly the same

  # Correlations of neighboring pixels, in DN^2
  #
  tCH = tCV = tCD = 0
  for if1 in range(1,num_files):
    for if2 in range(if1):
      temp_box = box2[if1,:,:] - box2[if2,:,:]

      # Run through twice if we have noise, otherwise once
      nrun = 2 if noise_corr else 1
      for icorr in range (nrun):
        # clipping
        cmin = pyIRC_percentile(temp_box,corr_mask,100*epsilon)
        cmax = pyIRC_percentile(temp_box,corr_mask,100*(1-epsilon))
        this_mask = numpy.where(numpy.logical_and(temp_box>cmin,temp_box<cmax),1,0) * corr_mask
        if numpy.sum(this_mask)<1: return [] # FAIL!
        # mean subtraction
        mean_of_temp_box = numpy.sum(temp_box*this_mask)/numpy.sum(this_mask)
        if subtr_corr and newMeanSubMethod: temp_box -= mean_of_temp_box

        # Correlations in horizontal and vertical directions
        maskCV = numpy.sum(this_mask[:-1,:]*this_mask[1:,:])
        maskCH = numpy.sum(this_mask[:,:-1]*this_mask[:,1:])
        CV = numpy.sum(this_mask[:-1,:]*this_mask[1:,:]*temp_box[:-1,:]*temp_box[1:,:])
        CH = numpy.sum(this_mask[:,:-1]*this_mask[:,1:]*temp_box[:,:-1]*temp_box[:,1:])
        if maskCH<1 or maskCV<1: return []
        CH /= maskCH
        CV /= maskCV

        # diagonal directions
        if not full_corr:
          maskCD1 = numpy.sum(this_mask[:-1,:-1]*this_mask[1:,1:])
          maskCD2 = numpy.sum(this_mask[:-1,1:]*this_mask[1:,:-1])
          CD1 = numpy.sum(this_mask[:-1,:-1]*this_mask[1:,1:]*temp_box[:-1,:-1]*temp_box[1:,1:])
          CD2 = numpy.sum(this_mask[:-1,1:]*this_mask[1:,:-1]*temp_box[:-1,1:]*temp_box[1:,:-1])
          if maskCD1<1 or maskCD2<1: return []
          CD1 /= maskCD1
          CD2 /= maskCD2
          CD = (CD1+CD2)/2.

        if leadtrailSub:
          maskCVx1 = numpy.sum(this_mask[:-1,:-4]*this_mask[1:,4:])
          maskCHx1 = numpy.sum(this_mask[:,:-5]*this_mask[:,5:])
          CVx1 = numpy.sum(this_mask[:-1,:-4]*this_mask[1:,4:]*temp_box[:-1,:-4]*temp_box[1:,4:])
          CHx1 = numpy.sum(this_mask[:,:-5]*this_mask[:,5:]*temp_box[:,:-5]*temp_box[:,5:])
          if maskCHx1<1 or maskCVx1<1: return []
          CHx1 /= maskCHx1
          CVx1 /= maskCVx1
          maskCVx2 = numpy.sum(this_mask[:-1,4:]*this_mask[1:,:-4])
          maskCHx2 = numpy.sum(this_mask[:,:-3]*this_mask[:,3:])
          CVx2 = numpy.sum(this_mask[:-1,4:]*this_mask[1:,:-4]*temp_box[:-1,4:]*temp_box[1:,:-4])
          CHx2 = numpy.sum(this_mask[:,:-3]*this_mask[:,3:]*temp_box[:,:-3]*temp_box[:,3:])
          if maskCHx2<1 or maskCVx2<1: return []
          CHx2 /= maskCHx2
          CVx2 /= maskCVx2
          CH -= (CHx1+CHx2)/2.
          CV -= (CVx1+CVx2)/2.
          #
          # correction of the diagonal directions
          if not full_corr:
            maskCDx1 = numpy.sum(this_mask[:-1,:-5]*this_mask[1:,5:])
            maskCDx2 = numpy.sum(this_mask[:-1,:-3]*this_mask[1:,3:])
            maskCDx3 = numpy.sum(this_mask[1:,:-5]*this_mask[:-1,5:])
            maskCDx4 = numpy.sum(this_mask[1:,:-3]*this_mask[:-1,3:])
            CDx1 = numpy.sum(this_mask[:-1,:-5]*this_mask[1:,5:]*temp_box[:-1,:-5]*temp_box[1:,5:])
            CDx2 = numpy.sum(this_mask[:-1,:-3]*this_mask[1:,3:]*temp_box[:-1,:-3]*temp_box[1:,3:])
            CDx3 = numpy.sum(this_mask[1:,:-5]*this_mask[:-1,5:]*temp_box[1:,:-5]*temp_box[1:,5:])
            CDx4 = numpy.sum(this_mask[1:,:-3]*this_mask[:-1,3:]*temp_box[1:,:-3]*temp_box[1:,3:])
            if maskCDx1<1 or maskCDx2<1 or maskCDx3<1 or maskCDx4<1: return []
            CDx1 /= maskCDx1
            CDx2 /= maskCDx2
            CDx3 /= maskCDx3
            CDx4 /= maskCDx4
            CD -= (CDx1+CDx2+CDx3+CDx4)/4.

        if subtr_corr and not newMeanSubMethod and not leadtrailSub:
          CH -= mean_of_temp_box**2
          CV -= mean_of_temp_box**2
        tCH += CH * (1 if icorr==0 else -1)
        tCV += CV * (1 if icorr==0 else -1)
        if not full_corr:
          if subtr_corr and not newMeanSubMethod and not leadtrailSub: CD -= mean_of_temp_box**2
          tCD += CD * (1 if icorr==0 else -1)

        if verbose:
          print ('pos =', if1, if2, 'iteration', icorr, 'cmin,cmax =', cmin, cmax)
          print ('Mask size', numpy.sum(this_mask), 'correlations =', maskCH, maskCV, 'data:', CH, CV)

        temp_box = box2Noise[if1,:,:] - box2Noise[if2,:,:]
        # end nested for loop
  #
  # Normalize covariances. Note that taking the difference of 2 frames doubled the covariance
  # matrix, so we have introduced cov_clip_corr
  xi = scipy.stats.norm.ppf(1-epsilon)
  cov_clip_corr = (1. - numpy.sqrt(2./numpy.pi)*xi*numpy.exp(-xi*xi/2.)/(1.-2.*epsilon) )**2
  tCH /= num_files*(num_files-1)*cov_clip_corr
  tCV /= num_files*(num_files-1)*cov_clip_corr
  if not full_corr: tCD /= num_files*(num_files-1)*cov_clip_corr

  # if we don't need full correlations, exit now
  if not full_corr:
    return [numpy.sum(this_mask), med2, var2, tCH, tCV, tCD]

  # Curvature information (for 2nd order NL coefficient)
  if (tslices[-1]!=tslices[-2]):
    if subtr_href:
      for f in range(num_files):
        box1[f,:,:] += lightref[f,nt+1]
    boxD = region_cube[0:num_files,-2,:,:] - region_cube[0:num_files,-1,:,:]\
           - (tslices[-1]-tslices[-2])/float(tslices[1]-tslices[0])*box1
           # difference map
    if subtr_href:
      for f in range(num_files):
        box1[f,:,:] -= lightref[f,nt+1]
        boxD[f,:,:] -= (tslices[-1]-tslices[-2])/float(tslices[1]-tslices[0]) * lightref[f,2*nt]
    fac0 = fac1 = 0
    for if1 in range(num_files):
      box1R = box1[if1,:,:]
      boxDR = boxD[if1,:,:]
      c1min = pyIRC_percentile(box1R, corr_mask, 100*epsilon)
      if c1min<=.5: c1min = .5   # should have no effect if successful, but prevents division by 0 if failure
      c1max = pyIRC_percentile(box1R, corr_mask, 100*(1-epsilon))
      cDmin = pyIRC_percentile(boxDR, corr_mask, 100*epsilon)
      cDmax = pyIRC_percentile(boxDR, corr_mask, 100*(1-epsilon))
      this_file_mask = numpy.where(numpy.logical_and(box1R>c1min, numpy.logical_and(box1R<c1max,
        numpy.logical_and(boxDR>cDmin, boxDR<cDmax))), corr_mask, 0)
      fac0 += numpy.sum(this_file_mask*boxDR)
      fac1 += numpy.sum(this_file_mask*box1R)
    if fac1<.5: return [] # FAIL!
    frac_dslope = fac0/fac1/((tslices[-1]-tslices[-2])/float(tslices[1]-tslices[0]))
  else:
    frac_dslope = 0.
  if verbose: print ('frac_dslope =', frac_dslope)

  if verbose:
    print ('Group 1 ->', med1, var1)
    print ('Group 2 ->', med2, var2)
    print ('correlations in Group 2:', tCH, tCV)
    print ('factors used: xi =', xi, ', cov_clip_corr =', cov_clip_corr)

  # Get alpha-corrected gains
  out = gain_alphacorr(gain_raw, tCH, tCV, med2)
  if tslices[1]>=tslices[-1] and len(out)<1:
    return [numpy.sum(this_mask), gain_raw, gain_raw, gain_raw, 0., 0., 0., med2/gain_raw/(tslices[1]-tslices[0]), 0., tCH, tCV]
  if len(out)<1: return [] # FAIL!
  gain_acorr = out[0]
  aH = out[1]
  aV = out[2]

  if tslices[1]>=tslices[-1]:
    return [numpy.sum(this_mask), gain_raw, gain_acorr, gain_acorr, aH, aV, 0., med2/gain_acorr/(tslices[1]-tslices[0]), 0., tCH, tCV]

  out = gain_alphabetacorr(gain_raw, tCH, tCV, med2, frac_dslope, [t-treset for t in tslices])
  if len(out)<1: return [] # FAIL!
  gain_abcorr = out[0]
  aH = out[1]
  aV = out[2]
  beta = out[3]
  I = out[4]

  return [numpy.sum(this_mask), gain_raw, gain_acorr, gain_abcorr, aH, aV, beta, I, 0., tCH, tCV]

# Under construction correlation functions for charge diffusion measurements
# many parts drawn from basic, might want option to do the usual 3x3 vs 5x5 but
# this could be done in principle just with the usual basic function?
# There's probably a better way of writing this...?!
def corr_5x5(region_cube, dark_cube, tslices, lightref, darkref, ctrl_pars, verbose):

  # Settings:
  newMeanSubMethod = True     # use False only for test/debug
  leadtrailSub = True         # subtract leading & trailing (by +/-4 pix) from horiz & vert correlations

  g_ptile = 75.               # percentile use for inter-quantile range for variance (default: 75, giving standard IQR)

  # Extract basic parameters
  num_files = region_cube.shape[0]-1
  nt = region_cube.shape[1]
  dy = region_cube.shape[2]
  dx = region_cube.shape[3]
  npix = dx*dy
  if nt!=len(tslices):
    print ('Error in pyirc.corr_5x5: incompatible number of time slices')
    exit()
  if verbose: print ('nfiles = ',num_files,', ntimes = ',nt,', dx,dy=',dx,dy)
  treset = 0
  if hasattr(ctrl_pars,'reset_frame'): treset = ctrl_pars.reset_frame

  # First get correlation parameters
  epsilon = .01
  if hasattr(ctrl_pars,'epsilon'): epsilon = ctrl_pars.epsilon
  subtr_corr = True
  if hasattr(ctrl_pars,'subtr_corr'): subtr_corr = ctrl_pars.subtr_corr
  noise_corr = True
  if hasattr(ctrl_pars,'noise_corr'): noise_corr = ctrl_pars.noise_corr
  if verbose: print ('corr pars =', epsilon, subtr_corr, noise_corr)
  #

  # Reference pixel subtraction?
  subtr_href = True
  if hasattr(ctrl_pars,'subtr_href'): subtr_href = ctrl_pars.subtr_href

  # lead-trail subtraction for IPC correlations?
  if hasattr(ctrl_pars,'leadtrailSub'): leadtrailSub = ctrl_pars.leadtrailSub

  # quantile for variance?
  if hasattr(ctrl_pars,'g_ptile'): g_ptile = ctrl_pars.g_ptile

  # Get means and variances at the early and last slices
  # (i.e. 1-point information)
  gauss_iqr_in_sigmas = scipy.stats.norm.ppf(g_ptile/100.)*2  # about 1.349 for g_ptile=75.
  box1 = region_cube[0:num_files,0,:,:] - region_cube[0:num_files,1,:,:]
  box2 = region_cube[0:num_files,0,:,:] - region_cube[0:num_files,-1,:,:]
  box2Noise = dark_cube[0:num_files,0,:,:] - dark_cube[0:num_files,-1,:,:]
  #
  if subtr_href:
    for f in range(num_files):
      if verbose: print ('lightref.shape=',lightref.shape, 'subtr ->', lightref[f,nt+1], lightref[f,2*nt-1], darkref[f,2*nt-1])
      box1[f,:,:] -= lightref[f,nt+1]
      box2[f,:,:] -= lightref[f,2*nt-1]
      box2Noise[f,:,:] -= darkref[f,2*nt-1]
  mean1 = numpy.mean(box1, axis=0)
  mean2 = numpy.mean(box2, axis=0)
  med1 = numpy.median(mean1)
  med2 = numpy.median(mean2)
  var1 = 0
  var2 = 0
  corr_mask = region_cube[-1,0,:,:]

  C_shift_mean = numpy.zeros((dy,dx))
  tC_all = numpy.zeros((dy,dx))

  for if1 in range(1,num_files):
    for if2 in range(if1):
      temp_box = box1[if1,:,:] - box1[if2,:,:]
      iqr1 = pyIRC_percentile(temp_box,corr_mask,g_ptile) - pyIRC_percentile(temp_box,corr_mask,100-g_ptile)
      temp_box = box2[if1,:,:] - box2[if2,:,:]
      iqr2 = pyIRC_percentile(temp_box,corr_mask,g_ptile) - pyIRC_percentile(temp_box,corr_mask,100-g_ptile)
      var1 += (iqr1/gauss_iqr_in_sigmas)**2/2.
      var2 += (iqr2/gauss_iqr_in_sigmas)**2/2.
      if verbose: print ('Inner loop,', if1, if2, temp_box.shape)

  var1 /= num_files*(num_files-1)/2.
  var2 /= num_files*(num_files-1)/2.
  if var2<=var1 and tslices[1]<tslices[-1]: return [] # FAIL!

  # Correlations of neighboring pixels, in DN^2
  #
  for if1 in range(1,num_files):
    for if2 in range(if1):
      temp_box = box2[if1,:,:] - box2[if2,:,:] 

      # Run through twice if we have noise, otherwise once
      nrun = 2 if noise_corr else 1
      if verbose: print("if1,if2=", if1, if2, " nrun: ",nrun)
      for icorr in range (nrun):
        # clipping
        cmin = pyIRC_percentile(temp_box,corr_mask,100*epsilon)
        cmax = pyIRC_percentile(temp_box,corr_mask,100*(1-epsilon))
        this_mask = numpy.where(numpy.logical_and(temp_box>cmin,temp_box<cmax),1,0) * corr_mask
        if numpy.sum(this_mask)<1: return [] # FAIL!
        # mean subtraction
        mean_of_temp_box = numpy.sum(temp_box*this_mask)/numpy.sum(this_mask)
        if subtr_corr and newMeanSubMethod: temp_box -= mean_of_temp_box

        # Correlations in all directions
        #masktmp = correlate2d(this_mask, this_mask,mode='same')
        #C_all = correlate2d(this_mask*temp_box, this_mask*temp_box, mode='same')
        dy2 = dy//2; dx2 = dx//2
        masktmp = fftconvolve(this_mask, numpy.flip(this_mask),mode='full')[dy2:-dy2+1,dx2:-dx2+1]
        C_all = fftconvolve(this_mask*temp_box, numpy.flip(this_mask*temp_box), mode='full')[dy2:-dy2+1,dx2:-dx2+1]

        if numpy.any(masktmp<1): return []

        C_all /= masktmp

        if leadtrailSub:
          C_pos_shift = numpy.zeros_like(C_all)
          C_neg_shift = numpy.zeros_like(C_all)

          C_pos_shift[:,:-8]=C_all[:,8:] #values of the correlation matrix 8 columns to the right
          C_neg_shift[:,8:]=C_all[:,:-8] #values of the correlation matrix 8 columns to the left

          """The 8 columns at the right edge just take the negative shift values, 
             the 8 columns at the left edge just take the positive shift values,
             and in the middle the mean of the two shifts is computed:
          """
          C_shift_mean[:, 8:-8] = numpy.mean([C_pos_shift[:, 8:-8], C_neg_shift[:, 8:-8]], axis=0)
          C_shift_mean[:, :8] = C_pos_shift[:, :8]
          C_shift_mean[:, -8:] = C_neg_shift[:, -8:]

          C_all = C_all - C_shift_mean

        #need to update the lines below to use C_all
        if subtr_corr and not newMeanSubMethod and not leadtrailSub:
          C_all -= mean_of_temp_box**2

        tC_all += C_all * (1 if icorr==0 else -1)

        if verbose:
          print ('pos =', if1, if2, 'iteration', icorr, 'cmin,cmax =', cmin, cmax)
          # Below needs to be adjusted
          #print ('Mask size', numpy.sum(this_mask), 'correlations =', maskCH, maskCV, 'data:', CH, CV)

        temp_box = box2Noise[if1,:,:] - box2Noise[if2,:,:]
        # end nested for loop

  #
  # Normalize covariances. Note that taking the difference of 2 frames doubled the covariance
  # matrix, so we have introduced cov_clip_corr
  xi = scipy.stats.norm.ppf(1-epsilon)
  cov_clip_corr = (1. - numpy.sqrt(2./numpy.pi)*xi*numpy.exp(-xi*xi/2.)/(1.-2.*epsilon) )**2
  tC_all /= num_files*(num_files-1)*cov_clip_corr

  # extract 5x5 matrix in the center of tC_all here:
  # hard-coded to return only 5x5 arrays, we should add option to specify
  # Find the "center" of this array
  if (dy%2==0):
    c_y=dy//2
  else:
    c_y=dy/2 - 1
  if (dx%2==0):
    c_x=dx//2
  else:
    c_x=dx/2 - 1
  tC_all_5x5 = tC_all[c_y-3:c_y+2,c_x-3:c_x+2]
  decenter_tC_all = decenter(tC_all_5x5)  # Might come in handy
  if verbose: print('tCH, tCV: ', decenter_tC_all[0,1], decenter_tC_all[1,0])

  # Return the correlations
  return [numpy.sum(this_mask), med2, var1, var2, tC_all_5x5]

# Routine to obtain statistical properties of a region of the detector across many time slices
#
# Inputs:
# lightfiles = list of light files
# darkfiles = list of dark files
# formatpars = format parameters
# box = list [xmin, xmax, ymin, ymax]
# tslices = list [tmin, tmax, deltas ...] (Python format -- xmax, ymax, tmax not included)
#   if no deltas (tslices length 2) then compute everything; if deltas specified, then only compute
#   correlations at the specified deltas (e.g. [1,3] for delta t = 1 or 3)
# sensitivity_spread_cut = for good pixels (typically 0.1)
# ctrl_pars = parameters for basic
#
# Each data[ti,tj,:] contains:
#   [number of good pixels, median, variance, tCH, tCV, tCD]
def corrstats(lightfiles, darkfiles, formatpars, box, tslices, sensitivity_spread_cut, ctrl_pars):

  # make copy of ctrl_pars, but force 5th element to be False
  ctrl_pars2 = copy.copy(ctrl_pars)
  ctrl_pars2.full_corr = False

  tmin = tslices[0]; tmax = tslices[1]; nt = tmax-tmin
  # build cube of good pixels, medians, variances, correlations
  data = numpy.zeros((nt,nt,6))
  # and get mask (last 'time' slice) -- only thing we are extracting from region_cube_X
  region_cube_X = pixel_data(lightfiles, formatpars, box[:4], [tmin,tmax-1,tmax-1,tmax-1], [sensitivity_spread_cut, True], False)

  # Get list of (good pix, median, var, cov_H, cov_V)
  for ti in range(nt-1):
    for tj in range(ti+1,nt):
      if tslices[2:]==[] or tj-ti in tslices[2:] or tj-ti==nt-1:
        t1 = tmin+ti
        t2 = tmin+tj
        tarray = [t1,t2,t2,t2]
        lightref = ref_array_block(lightfiles, formatpars, box[2:4], tarray, False)
        darkref = ref_array_block(darkfiles, formatpars, box[2:4], tarray, False)
        if not ctrl_pars.subtr_href:
          lightref[:,:] = 0.
          darkref[:,:] = 0.
        region_cube = pixel_data(lightfiles, formatpars, box[:4], tarray, [sensitivity_spread_cut, False], False)
        dark_cube = pixel_data(darkfiles, formatpars, box[:4], tarray, [sensitivity_spread_cut, False], False)
        # switch to the mask from above
        region_cube[-1,:,:,:] = region_cube_X[-1,:,:,:]
        dark_cube[-1,:,:,:] = region_cube_X[-1,:,:,:]
        B = basic(region_cube, dark_cube, tarray, lightref, darkref, ctrl_pars2, False)
        if len(B)==6: data[ti,tj,:] = numpy.asarray(B)
        # print (t1, t2, data[ti,tj,:], len(B))

  return data

# Routine to characterize of a region of the detector across many time slices
#
# Inputs:
# lightfiles = list of light files
# darkfiles = list of dark files
# formatpars = format parameters
# box = list [xmin, xmax, ymin, ymax]
# tslices = list [tmin, tmax, dt1, dt2] (Python format -- xmax, ymax, tmax not included)
#   correlations at the specified dt's are used (e.g. [1,5])
# sensitivity_spread_cut = for good pixels (typically 0.1)
# ctrl_pars = parameters for basic
# addInfo = additional information (sometimes needed)
# 
# return value is [isgood (1/0), g, aH, aV, beta, I, aD, da (residual)]
#
def polychar(lightfiles, darkfiles, formatpars, box, tslices, sensitivity_spread_cut, ctrl_pars, addInfo):

  # Check whether we have non-linearity information
  if ctrl_pars.use_allorder:
    if len(addInfo)<3:
      print ('Error: polychar: not enough fields in addInfo')
      return []

  # Check time range
  if len(tslices)<4:
    print ('Error: polychar: not enough data', tslices)
    return []
  if tslices[2]>=tslices[3] or tslices[3]>=tslices[1]-tslices[0] or tslices[1]-tslices[0]<3:
    print ('Error: polychar: invalid slices range', tslices)
    return []

  # Get correlation function data (including adjacent steps))
  data = corrstats(lightfiles, darkfiles, formatpars, box, tslices+[1], sensitivity_spread_cut, ctrl_pars)

  # check if this is good
  nt = tslices[1]-tslices[0]
  for ti in range(nt-1):
    for tj in range(ti+1,nt):
      if data[ti,tj,0]==0 and tj-ti in [1,tslices[2],tslices[3]]:
        return [0,0,0,0,0,0]

  # Determine whether we are applying corrections
  applyCorr = False
  if len(addInfo)>=2:
    applyCorr = True
    typeCorr = addInfo[0]
    ipnl = addInfo[1]
    sBFE = numpy.shape(ipnl)[0]//2

  # Fit of differences as a function of slice number
  # slope = -2*beta*I^2/g
  # intercept = (I - beta I^2)/g
  npts = tslices[1]-tslices[0]-1
  diff_frames = numpy.zeros((npts))
  for j in range(npts):
    diff_frames[j] = data[j,j+1,1] # median from frame tslices[0]+j -> tslices[0]+j+1
  slopemed, icpt = numpy.linalg.lstsq(numpy.vstack([numpy.array(range(npts)) + tslices[0]-ctrl_pars.reset_frame,
                   numpy.ones(npts)]).T, diff_frames, rcond=-1)[0]
  # If using 'allorder', let's subtract out the higher-order terms:
  if ctrl_pars.use_allorder:
    xr = numpy.array(range(npts)) + tslices[0]-ctrl_pars.reset_frame
    i=100; err=10;
    etarget = 1e-9*numpy.abs(icpt)
    while i>=0 and err>etarget:
      I__g = icpt - 0.5*slopemed
      diff_frames_reduced = diff_frames.copy()
      icpt_old = icpt
      slopemed_old = slopemed
      for j in range(3, swi.p+1):
        diff_frames_reduced -= addInfo[2][j-2]*((xr+1)**j - xr**j) * (icpt-slopemed*.5)**j
      slopemed, icpt = numpy.linalg.lstsq(numpy.vstack([xr, numpy.ones(npts)]).T, diff_frames_reduced, rcond=-1)[0]
      err = numpy.sqrt( (icpt-icpt_old)**2 + (slopemed-slopemed_old)**2 )
      if i==0:
        print ('higher order loop failed to converge {:12.5E} vs {:12.5E} (target)', err, etarget)
        return []

  # Difference of correlation functions
  #
  # Cdiff = I/g^2 * ((1-4a)^2 + 2aH^2 + 2aV^2) * t_{bd} - 4(1-8a)beta I^2/g^2 * (t_{ad}t_d - t_{ab}t_b + (e-1)/2*t_{bd})
  # where e = npts2 is number of bins averaged together
  #
  # and horizontal and vertical cross-correlations
  # CH = 2 I t_{ab} / g^2 * ( 1-4a - 4 beta (I t_b + 1/2 + (e-1)/2*I) ) * aH
  # CV = 2 I t_{ab} / g^2 * ( 1-4a - 4 beta (I t_b + 1/2 + (e-1)/2*I) ) * aV
  #
  npts2 = tslices[1]-tslices[0]-tslices[3]
  Cdiff = CV = CH = CD = 0.
  for j in range(npts2):
    Cdiff += data[j,j+tslices[3],2] - data[j,j+tslices[2],2]
    CH += data[j,j+tslices[3],3]
    CV += data[j,j+tslices[3],4]
    CD += data[j,j+tslices[3],5]
  Cdiff /= npts2; CH /= npts2; CV /= npts2; CD /= npts2

  # initialize with no IPC or NL
  alphaH = alphaV = alphaD = alpha = beta = 0.
  da = 1.
  # dummy initializations; these get over-written before they are used
  I = g = 1.
  Cdiffcorr = 0.
  iCycle = 0; nCycle=100
  while iCycle<nCycle:
    alphaH_old = alphaH; alphaV_old = alphaV; alphaD_old=alphaD; g_old=g # to track convergence

    # Get combination of I and gain from difference of correlation functions
    tbrack = tslices[3]*(tslices[0]+tslices[3]-ctrl_pars.reset_frame) - tslices[2]*(tslices[0]+tslices[2]-ctrl_pars.reset_frame)\
             + (npts2-1)/2.0*(tslices[3]-tslices[2])
    I__g2 = (Cdiff - Cdiffcorr + 4.*(1.-8.*alpha)*beta*I**2/g**2*tbrack) / (tslices[3]-tslices[2]) / ( (1.-4*alpha-4*alphaD)**2 + 2*alphaH**2+2*alphaV**2 + 4*alphaD**2 )

    # Now use slopemed = -2 beta I^2/g, icpt = (I - beta I^2)/g, and I/g^2 to solve for I, beta, and g
    g = (icpt - slopemed/2.)/I__g2
    I = I__g2 * g**2
    beta = -g*slopemed/2./I**2

    # Corrections to horiz. and vert. IPC
    #
    CHcorr = CVcorr = CDcorr = 0.
    if applyCorr:
      if typeCorr.lower() == 'bfe':
        CHcorr = (ipnl[sBFE,sBFE+1]+ipnl[sBFE,sBFE-1])/2. * (I/g*tslices[3])**2
        CVcorr = (ipnl[sBFE+1,sBFE]+ipnl[sBFE-1,sBFE])/2. * (I/g*tslices[3])**2
        CDcorr = (ipnl[sBFE+1,sBFE+1]+ipnl[sBFE+1,sBFE-1]+ipnl[sBFE-1,sBFE+1]+ipnl[sBFE-1,sBFE-1])/4. * (I/g*tslices[3])**2
        Cdiffcorr = ipnl[sBFE,sBFE] * (I/g)**2*(tslices[3]**2-tslices[2]**2)
      if typeCorr.lower() == 'nlipc':
        CHcorr = (ipnl[sBFE,sBFE+1]+ipnl[sBFE,sBFE-1])/2. * (I/g)**2*tslices[3]*(tslices[0]+tslices[3]+(npts2-1)*0.5)*2
        CVcorr = (ipnl[sBFE+1,sBFE]+ipnl[sBFE-1,sBFE])/2. * (I/g)**2*tslices[3]*(tslices[0]+tslices[3]+(npts2-1)*0.5)*2
        CDcorr = (ipnl[sBFE+1,sBFE+1]+ipnl[sBFE+1,sBFE-1]+ipnl[sBFE-1,sBFE+1]+ipnl[sBFE-1,sBFE-1])/4. * (I/g)**2*tslices[3]*(tslices[0]+tslices[3]+(npts2-1)*0.5)*2
        Cdiffcorr = ipnl[sBFE,sBFE] * (I/g)**2*( (tslices[0]+tslices[3])*tslices[3] - (tslices[0]+tslices[2])*tslices[2]
                      + (tslices[3]-tslices[2])*(npts2-1)*0.5)
      
      # apply corrections from ftsolve
      if ctrl_pars.fullnl and typeCorr.lower() == 'bfe':
        beta_cm = beta
        if ctrl_pars.use_allorder: beta_cm = -addInfo[2]/g**numpy.linspace(1,swi.p-1,num=swi.p-1)
        if Test_SubBeta: beta_cm = beta
        t0 = tslices[0]-ctrl_pars.reset_frame
        CF_BigStep = solve_corr_many(ipnl, 21, I, g, beta_cm, 0., [t0, t0+tslices[3], t0, t0+tslices[3], npts2],
          [alphaV, alphaH, alphaD], [0.,0.,0.], sBFE)
        CF_SmallStep = solve_corr_many(ipnl, 21, I, g, beta_cm, 0., [t0, t0+tslices[2], t0, t0+tslices[2], npts2],
          [alphaV, alphaH, alphaD], [0.,0.,0.], sBFE)
        Cdiffcorr = CF_BigStep[sBFE,sBFE] - CF_SmallStep[sBFE,sBFE] - (
          I/g**2*((1-4*alpha-4*alphaD)**2+2*alphaH**2+2*alphaV**2+4*alphaD**2)*(tslices[3]-tslices[2])
          -4*(1-8*alpha)*beta*I**2/g**2*tbrack)
        ad3 = tslices[0]+tslices[3]-ctrl_pars.reset_frame
        CHcorr = (CF_BigStep[sBFE,sBFE+1]+CF_BigStep[sBFE,sBFE-1])/2. - (
          2.*I/g**2*tslices[3]*(1.-4*alpha-4*alphaD-4*beta*(I*ad3+(npts2-1)*0.5*I+0.5))*alphaH + 4.*I/g**2*tslices[3]*alphaV*alphaD)
        CVcorr = (CF_BigStep[sBFE+1,sBFE]+CF_BigStep[sBFE-1,sBFE])/2. - (
          2.*I/g**2*tslices[3]*(1.-4*alpha-4*alphaD-4*beta*(I*ad3+(npts2-1)*0.5*I+0.5))*alphaV + 4.*I/g**2*tslices[3]*alphaH*alphaD)
        CDcorr = (CF_BigStep[sBFE+1,sBFE+1]+CF_BigStep[sBFE-1,sBFE+1]+CF_BigStep[sBFE+1,sBFE-1]+CF_BigStep[sBFE-1,sBFE-1])/4. - (
          2.*I/g**2*tslices[3]*(1.-4*alpha-4*alphaD)*alphaD + 2.*I/g**2*tslices[3]*alphaH*alphaV)

    factor = 2.*I__g2*tslices[3] * ( 1.-4.*alpha - 4.*alphaD - 4.*beta*( I*(tslices[0]+tslices[3]-ctrl_pars.reset_frame+(npts2-1.)/2.) +0.5) )
    factor_raw = 2.*I__g2*tslices[3]
    alphaH = (CH - CHcorr - 2.*alphaV*alphaD*factor_raw)/factor
    alphaV = (CV - CVcorr - 2.*alphaH*alphaD*factor_raw)/factor
    alphaD = ( (CD - CDcorr)/factor_raw - alphaH*alphaV) / (1.-4.*alpha-4.*alphaD)
    alpha = (alphaH+alphaV)/2.
    da = numpy.abs(alphaH_old-alphaH) + numpy.abs(alphaV_old-alphaV) + numpy.abs(alphaD_old-alphaD)
    dg = numpy.abs(g_old-g)
    iCycle+=1
    if iCycle<nCycle-2 and da<1e-8 and dg<1e-8: iCycle=nCycle-2 # fast exit from loop

  return [1, g, alphaH, alphaV, beta, I, alphaD, da]

# Routines to compute the BFE coefficients
#
# Inputs:
# region_cube = 4D array of the region of interest (order of indices: file, timeslice, y-ymin, x-xmin)
# tslices = list of time slices
# basicinfo = output from basic (incl. gains, IPC, NL)
# ctrl_pars_bfe = parameters to control BFE determination
#   ctrl_pars_bfe.epsilon = cut fraction (default to 0.01)
#   ctrl_pars_bfe.treset = reset frame (default to 0)
#   ctrl_pars_bfe.BSub = baseline subtraction? (default to True)
#   ctrl_pars_bfe.vis = has visible? (default to False)
#   ctrl_pars_bfe.Phi = omega*p2/(1+omega) kernel (only used if ctrl_pars_bfe.vis is true)
# verbose = True or False (recommend True only for debugging)
#
# output is a fsBFE x fsBFE (default: 5x5) BFE kernel in inverse electrons
#
def bfe(region_cube, tslices, basicinfo, ctrl_pars_bfe, verbose):
  N = 21 # <-- size for ftsolve

  # Extract parameters from basicinfo
  gain =   basicinfo[swi.g]
  aH =     basicinfo[swi.alphaH]
  aV =     basicinfo[swi.alphaV]
  beta =   basicinfo[swi.beta]
  I =      basicinfo[swi.I]

  # Extract basic parameters
  num_files = region_cube.shape[0]-1
  nt = region_cube.shape[1]
  dy = region_cube.shape[2]
  dx = region_cube.shape[3]
  npix = dx*dy
  if (nt!=len(tslices)):
    print ('Error in basic: incompatible number of time slices')
    exit()
  if verbose: print ('nfiles = ',num_files,', ntimes = ',nt,', dx,dy=',dx,dy)
  treset = 0
  if hasattr(ctrl_pars_bfe,'treset'): treset = ctrl_pars_bfe.treset

  # for visible flats
  hasvis = False
  if hasattr(ctrl_pars_bfe,'vis'):
    if ctrl_pars_bfe.vis:
      hasvis = True
      normPhi = numpy.sum(ctrl_pars_bfe.Phi) # this is omega/(1+omega)
      omega = normPhi / (1-normPhi)
      p2 = numpy.zeros_like(ctrl_pars_bfe.Phi)
      if numpy.abs(normPhi)>1e-49: p2 = ctrl_pars_bfe.Phi / normPhi # this prevents an exception if omega=0
      p2 = pad_to_N(p2,N) # still centered

  # BFE kernel size:
  # sBFE = range; fsBFE = full size
  sBFE = swi.s
  fsBFE = 2*sBFE+1
  sBFE_out = sBFE
  fsBFE_out = fsBFE

  # replace beta with a scalar value if necessary
  # note beta[0] is now 2nd order coef (in DN^-1) is to be converted to beta (in e^-1) and has opposite sign
  if ctrl_pars_bfe.fullnl and ctrl_pars_bfe.use_allorder: beta = -beta/gain**numpy.linspace(1,swi.p-1,num=swi.p-1)
  if ctrl_pars_bfe.fullnl and ctrl_pars_bfe.use_allorder and Test_SubBeta: beta = beta[0]

  # Baseline subtraction -- requires bigger box
  BSub = True
  if hasattr(ctrl_pars_bfe,'BSub'): BSub = ctrl_pars_bfe.BSub
  if BSub:
    sBFE = max(sBFE_out, 10)
    fsBFE = 2*sBFE+1
    pad = 5 # Number of pixels in corr. fcn. to take for the baseline on each side in each row

  # Cut fraction and correction
  epsilon = .01
  if hasattr(ctrl_pars_bfe,'epsilon'): epsilon = ctrl_pars_bfe.epsilon
  xi = scipy.stats.norm.ppf(1-epsilon)
  cov_clip_corr = (1. - numpy.sqrt(2./numpy.pi)*xi*numpy.exp(-xi*xi/2.)/(1.-2.*epsilon) )**2

  # Build the two slices to correlate
  box1 = region_cube[0:num_files,0,:,:] - region_cube[0:num_files,1,:,:]
  box3 = region_cube[0:num_files,-2,:,:] - region_cube[0:num_files,-1,:,:]
  corr_mask = region_cube[-1,0,:,:]

  # setup for BFE kernel
  numBFE = numpy.zeros((fsBFE,fsBFE))
  denBFE = numpy.zeros((fsBFE,fsBFE))

  # Loop over the flat pairs we are going to use
  for if1 in range(1,num_files):
    for if2 in range(if1):
      # Build slices and mask
      slice_ab = box1[if1,:,:] - box1[if2,:,:]
      slice_cd = box3[if1,:,:] - box3[if2,:,:]
      ab_min = pyIRC_percentile(slice_ab, corr_mask, 100*epsilon)
      ab_max = pyIRC_percentile(slice_ab, corr_mask, 100*(1-epsilon))
      cd_min = pyIRC_percentile(slice_cd, corr_mask, 100*epsilon)
      cd_max = pyIRC_percentile(slice_cd, corr_mask, 100*(1-epsilon))
      this_file_mask_ab = numpy.where(numpy.logical_and(slice_ab>ab_min, slice_ab<ab_max), corr_mask, 0)
      this_file_mask_cd = numpy.where(numpy.logical_and(slice_cd>cd_min, slice_cd<cd_max), corr_mask, 0)
      if verbose:
        print (if1, if2, slice_ab.shape, slice_cd.shape, numpy.sum(this_file_mask_ab), numpy.sum(this_file_mask_cd))

      # Mean subtraction
      slice_ab -= pyIRC_mean(slice_ab, this_file_mask_ab)
      slice_cd -= pyIRC_mean(slice_cd, this_file_mask_cd)
      # Set masked values to zero
      slice_ab *= this_file_mask_ab
      slice_cd *= this_file_mask_cd

      # Now get the correlation function ...
      # format is: numerator and denominator of C_{abcd}(2*sBFE-i,2*sBFE-j)
      for j in range(fsBFE):
        for i in range(fsBFE):
          abminX = 0
          abmaxX = dx
          abminY = 0
          abmaxY = dy
          if i>=sBFE:
            abmaxX += sBFE-i
          else:
            abminX += sBFE-i
          if j>=sBFE:
            abmaxY += sBFE-j
          else:
            abminY += sBFE-j
          cdminX = abminX + i - sBFE
          cdmaxX = abmaxX + i - sBFE
          cdminY = abminY + j - sBFE
          cdmaxY = abmaxY + j - sBFE

          # Add up contributions to the correlation function
          denBFE[j,i] += numpy.sum(this_file_mask_ab[abminY:abmaxY,abminX:abmaxX]*this_file_mask_cd[cdminY:cdmaxY,cdminX:cdmaxX])
          numBFE[j,i] += numpy.sum(slice_ab[abminY:abmaxY,abminX:abmaxX]*slice_cd[cdminY:cdmaxY,cdminX:cdmaxX])/2.
          # division by 2 since differencing two images doubles the answer

  BFEK = numBFE/(denBFE+1e-99)
  BFEK *= gain**2/(I**2*(tslices[1]-tslices[0])*(tslices[-1]-tslices[-2])*cov_clip_corr)

  # Baseline subtraction
  if BSub:
    for j in range(fsBFE):
      rowBL = ( numpy.mean(BFEK[j,0:pad]) + numpy.mean(BFEK[j,-pad:]) )/2.
      BFEK[j,:] -= rowBL

  # Implement cr_converge.
  if ctrl_pars_bfe.fullnl:
    avals = [basicinfo[swi.alphaV], basicinfo[swi.alphaH], basicinfo[swi.alphaD]]
    avals_nl = [0,0,0]
    sigma_a = 0
    tol = 1.e-11 #Pick a tolerance below which the two Crs are considered equal
    fsBFE_out = 2*sBFE_out+1
    observed_Cr = BFEK[sBFE-sBFE_out:sBFE+sBFE_out+1, sBFE-sBFE_out:sBFE+sBFE_out+1]
    BFEK_model = numpy.zeros((fsBFE_out,fsBFE_out))+1e-15
    element_diff = 10
    iters = 0
    while element_diff > tol and iters<=100:
        # Note: solve_corr takes centered things, decenters/calculates internally
        if hasvis:
          theory_Cr = solve_corr_vis(BFEK_model,N,I,gain,beta,sigma_a,[t-treset for t in tslices],avals,avals_nl,sBFE_out,omega,p2)\
            *((gain**2)/(I**2*(tslices[1]-tslices[0])*(tslices[-1]-tslices[-2])))
        else:
          theory_Cr = solve_corr(BFEK_model,N,I,gain,beta,sigma_a,[t-treset for t in tslices],avals,avals_nl)\
            *((gain**2)/(I**2*(tslices[1]-tslices[0])*(tslices[-1]-tslices[-2])))
        if numpy.isnan(theory_Cr).any():
            warnings.warn('BFE loop diverged, generated NaN')
            return numpy.zeros((fsBFE_out,fsBFE_out)) + numpy.nan
        difference = theory_Cr - observed_Cr
        element_diff = numpy.amax(abs(difference))
        BFEK_model -= difference[::-1,::-1]
        iters += 1
        if verbose: print(iter, BFEK_model)
        if iters>99:
           warnings.warn("WARNING: NL loop has iterated 100 times")
           return numpy.zeros((fsBFE_out,fsBFE_out)) + numpy.nan
    return BFEK_model

  else:
     # Corrections for classical non-linearity
     BFEK[sBFE,sBFE] += 2*(1-4*(aH+aV))*beta
     if sBFE>=1:
       BFEK[sBFE,sBFE+1] += 4*aH*beta
       BFEK[sBFE,sBFE-1] += 4*aH*beta
       BFEK[sBFE+1,sBFE] += 4*aV*beta
       BFEK[sBFE-1,sBFE] += 4*aV*beta
     return BFEK[sBFE-sBFE_out:sBFE+sBFE_out+1, sBFE-sBFE_out:sBFE+sBFE_out+1]

# Hot pixel identification
# Returns a tuple of hot pixels in the array that meet the following criteria:
# (*) apparent brightness in time slices up through tslices[-1] is assessed
# (*) in range from pars[0] .. pars[1] in last slice
# (*) repeatable to within a top-to-bottom error of pars[2] as a fraction of the
#       maximum signal (e.g. 0.1 for 10% repeatability)
# (*) isolation: if pars[3]>0, rejects pixels with neighbors that are at least pars[3] times
#       as bright as this pixel itself (e.g. 0.1 for 10% isolation)
#
def hotpix(darkfiles, formatpars, tslices, pars, verbose):

  # Build array for the dark cube
  ndarks = len(darkfiles)
  N = get_nside(formatpars)
  cube = numpy.zeros((ndarks,N,N))
  for f in range(ndarks):
    CDS = load_segment(darkfiles[f], formatpars, [0,N,0,N], [1,tslices[-1]], False)
    cube[f,:,:] = CDS[0,:,:] - CDS[1,:,:]

  # Extract information on the pixels
  this_hot = numpy.zeros((N,N))
  ave_cube = numpy.mean(cube, axis=0)
  d_cube = numpy.max(cube, axis=0) - numpy.min(cube, axis=0)
  if verbose:
    print ('time slices for hot pixel analysis ->', tslices)
    print (ave_cube)
    print ('->', ave_cube.shape)
    print (d_cube)
    print ('->', d_cube.shape)

  this_hot = numpy.where(numpy.logical_and(ave_cube>=pars[0], ave_cube<=pars[1]), 1, 0)

  # Isolation cut
  if verbose: print ('Start with', numpy.sum(this_hot), 'pixels before isolation cut')
  if pars[3]>0:
    C = 2
    M = numpy.ones((2*C+1,2*C+1))
    M[C,C]=0
    isolation_mask = scipy.ndimage.maximum_filter(ave_cube, footprint=M, mode='constant', cval=0)
    # Also avoid pixels that border on reference pixels
    this_hot[:4+C,:] = 0
    this_hot[-(4+C):,:] = 0
    this_hot[:,:4+C] = 0
    this_hot[:,-(4+C):] = 0
    this_hot *= numpy.where(isolation_mask<=pars[3]*ave_cube, 1, 0)

  if verbose: print ('Start with', numpy.sum(this_hot), 'pixels')
  for t in tslices[1:]:
    for f in range(ndarks):
      CDS = load_segment(darkfiles[f], formatpars, [0,N,0,N], [1,t], False)
      cube[f,:,:] = CDS[0,:,:] - CDS[1,:,:]
    d_cube = numpy.max(cube, axis=0) - numpy.min(cube, axis=0)
    this_hot *= numpy.where(d_cube<=pars[2]*ave_cube, 1, 0)
  if verbose: print (numpy.sum(this_hot))

  return numpy.where(this_hot>0)

# Return IPC data from a list of hot pixels.
# y, x = lists of hot pixel coordinates to use (probably selected from hotpix)
# darkfiles = list of dark files to use
# formatpars = format code for dark files
# tslices = list of time slices to report
# pars = parameters to control data selection
#        right now, if not empty:
#        pars[0] = numpy array map of non-linearity * gain values (units: DN^-1) (skip if not numpy.ndarray)
#        pars[1] = reference non-linearity to median stack of initial image? (T/F)
# verbose = T/F
#
# Returns data cube with three indices.
# data[jpix,jt,jpos] = signal (in DN) at position jpos (relative to hot pixel)
#                      at time slice jt
#                      of hot pixel jpix
#
# positions are: jpos=0 (center), 1 (right), 2 (up), 3 (left), 4 (down), 5-8 (diag, quadrants I-IV), 9 (bkgnd 5x5-3x3)
def hotpix_ipc(y, x, darkfiles, formatpars, tslices, pars, verbose):

  # Build array for the dark cube
  ndarks = len(darkfiles)
  N = get_nside(formatpars)
  cube = numpy.zeros((ndarks,N,N))

  nt = len(tslices)
  npix = len(x)
  data = numpy.zeros((npix,nt,10))

  # offset table
  dx = [0, 1, 0, -1, 0, 1, -1, -1, 1]
  dy = [0, 0, 1, 0, -1, 1, 1, -1, -1]

  # Perform nonlinearity correction?
  do_nonlin = False
  if len(pars)>=1:
    if type(pars[0]) is numpy.ndarray:
      do_nonlin = True
      m = pars[0]
      beta_gain = numpy.zeros((N,N))
      (ny1,nx1) = numpy.shape(m)
      kx1 = N//nx1; ky1 = N//ny1
      for i in range(nx1):
        for j in range(ny1):
          beta_gain[ky1*j:ky1*(j+1),kx1*i:kx1*(i+1)] = m[j,i]
      # now beta_gain is an NxN map of beta*gain
      if verbose: print ('beta*gain =', beta_gain)
  # baseline for NL correction is median image?
  medbaseline_nonlin = False
  if len(pars)>=2:
    medbaseline_nonlin = pars[1]

  # background mask
  bkmask = numpy.ones((5,5))
  bkmask[1:4,1:4]=0.
  fourmask = False
  if fourmask:
    bkmask[:,:]=0.
    bkmask[2,0] = bkmask[2,4] = bkmask[0,2] = bkmask[4,2] = 1.
  if verbose: print ('bkmask =', bkmask)
  # 16 ones and 9 zeros

  # now make data cube
  for jt in range(nt):
    for f in range(ndarks):
      CDS = load_segment(darkfiles[f], formatpars, [0,N,0,N], [1,tslices[jt]], False)
      cube[f,:,:] = CDS[0,:,:] - CDS[1,:,:]
      if do_nonlin:
        # non-linearity correction, if turned on
        cube_corr = cube[f,:,:]
        if medbaseline_nonlin:
          cube_corr = 2*scipy.ndimage.median_filter(CDS[0,:,:],size=3) - CDS[0,:,:] - CDS[1,:,:]
        cube[f,:,:] = cube[f,:,:]*(1.+beta_gain*cube_corr)
    medframe = numpy.median(cube, axis=0)
    if verbose: print ('med', numpy.shape(medframe), jt)
    for jpix in range(npix):
      for jpos in range(9):
        x_ = x[jpix] + dx[jpos]
        y_ = y[jpix] + dy[jpos]
        data[jpix,jt,jpos] = medframe[y_,x_]
      data[jpix,jt,9] = 25./16.*numpy.mean(bkmask*medframe[y[jpix]-2:y[jpix]+3, x[jpix]-2:x[jpix]+3])
      if fourmask: data[jpix,jt,9] *= 16./4.

  return data

# sliding median function
#
# takes in vectors x and y (length N) and percentile p.
# returns slope m such that p% of the data are below the line y = m*x.
#
# works by bisection in the 'guess' range mrange, with niter=64 iterations as default.
# default mrange is [-1,1] (appropriate for IPC uses)
#
# The pivot is not used but is here for forward compatibility; right now it assumes
# that the pivot point of the distribution is at x>0 (hence ordering in the bisection).
# In a future release if we need to change this the functionality is there.
#
def slidemed_percentile(x,y,p,mrange=[-1,1],niter=64,pivot='pos'):
  m1 = mrange[0]
  m2 = mrange[1]

  for k in range(niter):
    m = (m1+m2)/2.
    if numpy.nanpercentile(y-m*x,p)>0:
      m1=m
    else:
      m2=m
  return m

# Generates min and max range for a color bar
# based on inter-quartile range
def get_vmin_vmax(mydata, qext):
  Q1 = numpy.nanpercentile(mydata,25)
  Q2 = numpy.nanpercentile(mydata,75)
  return Q1-(Q2-Q1)*qext, Q2+(Q2-Q1)*qext
