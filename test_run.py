import sys
import time
import re
import numpy
import pyirc
import matplotlib
import matplotlib.pyplot as plt

class EmptyClass:
  pass

outstem = 'default_output'
use_cmap = 'gnuplot'

mydet = ''
lightfiles = []
darkfiles = []
formatpars = 1
nx = 32
ny = 32
tslices = [3,11,13,21]
tslicesM2a = []
tslicesM2b = []
tslicesM3 = []
fullref = True
sensitivity_spread_cut = .1
critfrac = 0.75
mychar = 'Basic'
hotpix = False
ref_for_hotpix_is_autocorr = False
hotpix_logtspace = False
hotpix_slidemed = False

# order parameters
s_bfe = 2     # order of BFE parameters
p_order = 0   # non-linearity polynomial table coefficients (table at end goes through order p_order)
              # set to zero to turn this off

# Parameters for basic characterization
basicpar = EmptyClass()
basicpar.epsilon = .01
basicpar.subtr_corr = True
basicpar.noise_corr = True
basicpar.reset_frame = 1
basicpar.subtr_href = True
basicpar.full_corr = True
basicpar.leadtrailSub = False
basicpar.g_ptile = 75.
basicpar.fullnl = False
basicpar.use_allorder = False

# Parameters for BFE
bfepar = EmptyClass()
bfepar.epsilon = .01
bfepar.treset = basicpar.reset_frame
bfepar.blsub = True
bfepar.fullnl = False

# Plotting parameters
narrowfig = False

# Read in information
config_file = sys.argv[1]
with open(config_file) as myf: content = myf.read().splitlines()
is_in_light = is_in_dark = False
maskX = [] # list of regions to mask
maskY = []
for line in content:
  # Cancellations
  m = re.search(r'^[A-Z]+\:', line)
  if m: is_in_light = is_in_dark = False

  # Searches for files -- must be first given the structure of this script!
  if is_in_light:
    m = re.search(r'^\s*(\S.*)$', line)
    if m: lightfiles += [m.group(1)]
  if is_in_dark:
    m = re.search(r'^\s*(\S.*)$', line)
    if m: darkfiles += [m.group(1)]

  # -- Keywords go below here --

  # Search for outputs
  m = re.search(r'^OUTPUT\:\s*(\S*)', line)
  if m: outstem = m.group(1)
  # Search for input files
  m = re.search(r'^LIGHT\:', line)
  if m: is_in_light = True
  m = re.search(r'^DARK\:', line)
  if m: is_in_dark = True
  # Format
  m = re.search(r'^FORMAT:\s*(\d+)', line)
  if m: formatpars = int(m.group(1))

  # Bin sizes
  m = re.search(r'^NBIN:\s*(\d+)\s+(\d+)', line)
  if m:
    nx = int(m.group(1))
    ny = int(m.group(2))

  # Characterization type (Basic or Advanced)
  m = re.search(r'^CHAR:\s*(\S+)', line)
  if m:
     mychar = m.group(1)
     if mychar.lower()=='advanced':
       m = re.search(r'^CHAR:\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)', line)
       if m:
         tchar1 = int(m.group(2))
         tchar2 = int(m.group(3))
         ncycle = int(m.group(4))
         ipnltype = m.group(5)
       else:
         print ('Error: insufficient arguments: ' + line + '\n')
         exit()

  # Time slices
  m = re.search(r'^TIME:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
  if m: tslices = [ int(m.group(x)) for x in range(1,5)]
  m = re.search(r'^TIME2A:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
  if m: tslicesM2a = [ int(m.group(x)) for x in range(1,5)]
  m = re.search(r'^TIME2B:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
  if m: tslicesM2b = [ int(m.group(x)) for x in range(1,5)]
  m = re.search(r'^TIME3:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
  if m: tslicesM3 = [ int(m.group(x)) for x in range(1,5)]
  #
  # reference time slice
  m = re.search(r'^TIMEREF:\s*(\d+)', line)
  if m: bfepar.treset = basicpar.reset_frame = int(m.group(1))

  # reference pixel subtraction
  m = re.search(r'^REF\s+OFF', line)
  if m: fullref = False

  # sensitivity spread cut
  m = re.search(r'^SPREAD:\s*(\S+)', line)
  if m: sensitivity_spread_cut = float(m.group(1))

  # variance parameters
  m = re.search(r'^QUANTILE:\s*(\S+)', line)
  if m: basicpar.g_ptile = float(m.group(1))
  # correlation parameters
  m = re.search(r'^EPSILON:\s*(\S+)', line)
  if m: bfepar.epsilon = basicpar.epsilon = float(m.group(1))
  m = re.search(r'^IPCSUB:\s*(\S+)', line)
  if m: basicpar.leadtrailSub = m.group(1).lower() in ['true', 'yes']

  # Other parameters
  m = re.search(r'^DETECTOR:\s*(\S+)', line)
  if m: mydet = m.group(1)
  m = re.search(r'^COLOR:\s*(\S+)', line)
  if m: use_cmap = m.group(1)

  # Classical non-linearity
  m = re.search(r'^NLPOLY:\s*(\S+)\s+(\S+)\s+(\S+)', line)
  if m:
    p_order = int(m.group(1))
    nlfit_ts = int(m.group(2))
    nlfit_te = int(m.group(3))

  m = re.search(r'^FULLNL:\s*(\S+)\s+(\S+)\s+(\S+)', line)
  if m:
    basicpar.fullnl = m.group(1).lower() in ['true', 'yes']
    bfepar.fullnl = m.group(2).lower() in ['true', 'yes']
    basicpar.use_allorder = m.group(3).lower() in ['true', 'yes']

  # Hot pixels
  # (adu min, adu max, cut stability, cut isolation)
  m = re.search(r'^HOTPIX:\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
  if m:
    hotpix = True
    hotpix_ADU_range = [ float(m.group(x)) for x in range(1,5)]
  #
  # change reference for hot pixels from last point to autocorr
  m = re.search(r'^HOTREF\s+AUTOCORR', line)
  if m: ref_for_hotpix_is_autocorr = True
  # log spacing for times?
  m = re.search(r'^HOTPIX\s+LOGTSPACE', line)
  if m: hotpix_logtspace = True
  # sliding median alpha method?
  m = re.search(r'^HOTPIX\s+SLIDEMED', line)
  if m: hotpix_slidemed = True

  # Mask regions by hand
  m = re.search(r'^MASK:\s*(\d+)\s+(\d+)', line)
  if m:
    maskX = maskX + [int(m.group(1))]
    maskY = maskY + [int(m.group(2))]

  # Control figures
  m = re.search(r'^NARROWFIG', line)
  if m: narrowfig = True

# set up array size parameters
pyirc.swi.addbfe(s_bfe)
pyirc.swi.addhnl(p_order)
print ('Number of output field per superpixel =', pyirc.swi.N)

# Check number of slices available
NTMAX = 16384
for f in lightfiles+darkfiles:
  nt = pyirc.get_num_slices(formatpars, f)
  if nt<NTMAX: NTMAX=nt

print ('Output will be directed to {:s}*'.format(outstem))
print ('Light files:', lightfiles)
print ('Dark files:', darkfiles)
print ('Time slices:', tslices, 'max=',NTMAX)
print ('Mask regions:', maskX, maskY)
# 
if len(lightfiles)!=len(darkfiles) or len(lightfiles)<2:
  print ('Failed: {:d} light files and {:d} dark files'.format(len(lightfiles), len(darkfiles)))
  exit()

# Additional parameters
# Size of a block
N = pyirc.get_nside(formatpars)
# Side lengths
dx = N//nx
dy = N//ny
# Pixels in a block
npix = dx*dy

# Make table of reference pixel corrections for Method 1
if fullref:
  lightref = pyirc.ref_array(lightfiles, formatpars, ny, tslices, False)
  darkref = pyirc.ref_array(lightfiles, formatpars, ny, tslices, False)
else:
  lightref = numpy.zeros((len(lightfiles), ny, 2*len(tslices)+1))
  darkref = numpy.zeros((len(darkfiles), ny, 2*len(tslices)+1))
basicpar.subtr_href = fullref

# more allocations
my_dim = pyirc.swi.N
full_info = numpy.zeros((ny,nx,my_dim))
is_good = numpy.zeros((ny,nx))

if p_order>0:
  # now coefficients for the info table
  # note that in 'abs' mode, the full_info[:,:,0] grid is not actually used, it
  #   is just there for consistency of the format
  # I moved this up here since we want to have these coefficients before the main program runs
  nlcubeX, nlfitX, nlderX, pcoefX = pyirc.gen_nl_cube(lightfiles, formatpars, [basicpar.reset_frame, nlfit_ts, nlfit_te], [ny,nx],
    full_info[:,:,0], 'abs', False)
  # fill in
  for iy in range(ny):
    for ix in range(nx):
      if pcoefX[1,iy,ix]!=0:
        full_info[iy,ix,pyirc.swi.Nbb] = -pcoefX[0,iy,ix]/pcoefX[1,iy,ix]
        for o in range(2,pyirc.swi.p+1):
          full_info[iy,ix,pyirc.swi.Nbb+o-1] = pcoefX[o,iy,ix]/pcoefX[1,iy,ix]**o
      else:
        full_info[iy,ix,pyirc.swi.Nbb] = -1e49 # error code

# Detector characterization data in a cube (basic characterization + BFE Method 1)
# Stdout calls are a progress indicator
#
print ('Method 1, progress of calculation:')
sys.stdout.write('|')
for iy in range(ny): sys.stdout.write(' ')
print ('| <- 100%')
sys.stdout.write('|')
for iy in range(ny):
  sys.stdout.write('*'); sys.stdout.flush()
  for ix in range(nx):
    region_cube = pyirc.pixel_data(lightfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], tslices,
                  [sensitivity_spread_cut, True], False)
    dark_cube = pyirc.pixel_data(darkfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], tslices,
                  [sensitivity_spread_cut, False], False)
    info = pyirc.basic(region_cube, dark_cube, tslices, lightref[:,iy,:], darkref[:,iy,:], basicpar, False)
    if len(info)>0:
      if mychar.lower()=='advanced':
        for iCycle in range(ncycle):
          bfeCoefs = pyirc.bfe(region_cube, tslices, info, bfepar, False)
          Cdata = pyirc.polychar(lightfiles, darkfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)],
                 [tslices[0], tslices[-1]+1, tchar1, tchar2], sensitivity_spread_cut, basicpar, [ipnltype, bfeCoefs]) # 1,3 or 9,19
          info[pyirc.swi.ind1:pyirc.swi.ind2] = numpy.asarray(Cdata[pyirc.swi.indp1:pyirc.swi.indp2])
      bfeCoefs = pyirc.bfe(region_cube, tslices, info, bfepar, False)
      info += bfeCoefs[0:2*pyirc.swi.s+1,0:2*pyirc.swi.s+1].flatten().tolist()
    else:
      info = numpy.zeros((pyirc.swi.Nbb)).tolist()

    if len(info)==pyirc.swi.Nbb:
      full_info[iy,ix,0:pyirc.swi.Nbb] = numpy.array(info)
    if info[0]>=npix*critfrac:
      is_good[iy,ix] = 1
    else:
      full_info[iy,ix,1:] = 0 # wipe out this super-pixel

print ('|')

# Mask regions
for mask_index in range(len(maskX)):
  ix = maskX[mask_index]
  iy = maskY[mask_index]
  is_good[iy,ix] = 0
  full_info[iy,ix,:] = 0 # wipe out this super-pixel

print (full_info.shape)
print ('Number of good regions =', numpy.sum(is_good))
mean_full_info = numpy.mean(numpy.mean(full_info, axis=0), axis=0)/numpy.mean(is_good)
print ('Mean info from good regions =', mean_full_info)
print ('')

# Non-linearity cube
ntSub = tslices[-1]
nlcube, nlfit, nlder = pyirc.gen_nl_cube(lightfiles, formatpars, ntSub, [ny,nx],
  full_info[:,:,pyirc.swi.beta]*full_info[:,:,pyirc.swi.I], 'dev', False)
thisOut = open(outstem+'_nl.txt', 'w')
for iy in range(ny):
  for ix in range(nx):
    thisOut.write('{:3d} {:3d} {:1d} {:9.6f} {:9.6f}'.format(iy,ix,int(is_good[iy,ix]),
      full_info[iy,ix,pyirc.swi.beta]*full_info[iy,ix,pyirc.swi.g]*1e6, full_info[iy,ix,pyirc.swi.g]))
    for it in range(ntSub):
      thisOut.write(' {:7.1f}'.format(nlcube[it,iy,ix]))
    thisOut.write('\n')
thisOut.close()
nlMean = numpy.zeros((ntSub))
nlFit = numpy.zeros((ntSub))
nlMeanDer = numpy.zeros((ntSub))
for it in range(ntSub):
  nlMean[it] = numpy.sum(is_good*nlcube[it,:,:])/numpy.sum(is_good)
  nlFit[it] = numpy.sum(is_good*nlfit[it,:,:])/numpy.sum(is_good)
  nlMeanDer[it] = numpy.sum(is_good*nlder[it,:,:])/numpy.sum(is_good)

# Multi-panel figure showing basic characterization
ar = nx/(ny+0.0)
spr = 2.2
matplotlib.rcParams.update({'font.size': 8})
F = plt.figure(figsize=(7,9))
S = F.add_subplot(3,2,1)
S.set_title(r'Good pixel map (%)')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(full_info[:,:,0]*100/(dx*dy), cmap=use_cmap, aspect=ar, interpolation='nearest', origin='lower',
  vmin=100*critfrac, vmax=100)
F.colorbar(im, orientation='vertical')
S = F.add_subplot(3,2,2)
S.set_title(r'Gain map $g$ (e/DN)')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
svmin, svmax = pyirc.get_vmin_vmax(full_info[:,:,pyirc.swi.g], spr)
im = S.imshow(full_info[:,:,pyirc.swi.g], cmap=use_cmap, aspect=ar, interpolation='nearest', origin='lower',
  vmin=svmin, vmax=svmax)
F.colorbar(im, orientation='vertical')
S = F.add_subplot(3,2,3)
S.set_title(r'IPC map $\alpha$ (%)')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
svmin, svmax = pyirc.get_vmin_vmax((full_info[:,:,pyirc.swi.alphaH]+full_info[:,:,pyirc.swi.alphaV])/2.*100., spr)
im = S.imshow((full_info[:,:,pyirc.swi.alphaH]+full_info[:,:,pyirc.swi.alphaV])/2.*100., cmap=use_cmap, aspect=ar,
  interpolation='nearest', origin='lower', vmin=svmin, vmax=svmax)
F.colorbar(im, orientation='vertical')
S = F.add_subplot(3,2,4)
S.set_title(r'Non-linearity map $\beta$ (ppm/e)')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
svmin, svmax = pyirc.get_vmin_vmax(full_info[:,:,pyirc.swi.beta]*1e6, spr)
im = S.imshow(full_info[:,:,pyirc.swi.beta]*1e6, cmap=use_cmap, aspect=ar, interpolation='nearest', origin='lower',
  vmin=svmin, vmax=svmax)
F.colorbar(im, orientation='vertical')
S = F.add_subplot(3,2,5)
S.set_title(r'Charge $It_{n,n+1}$ (e):')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
svmin, svmax = pyirc.get_vmin_vmax(full_info[:,:,pyirc.swi.I], spr)
im = S.imshow(full_info[:,:,pyirc.swi.I], cmap=use_cmap, aspect=ar, interpolation='nearest', origin='lower',
  vmin=svmin, vmax=svmax)
F.colorbar(im, orientation='vertical')
S = F.add_subplot(3,2,6)
S.set_title(r'IPNL $[K^2a+KK^\prime]_{0,0}$ (ppm/e):')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
svmin, svmax = pyirc.get_vmin_vmax(full_info[:,:,pyirc.swi.ker0]*1e6, spr)
im = S.imshow(full_info[:,:,pyirc.swi.ker0]*1e6, cmap=use_cmap, aspect=ar, interpolation='nearest', origin='lower',
  vmin=svmin, vmax=svmax)
F.colorbar(im, orientation='vertical')
F.set_tight_layout(True)
F.savefig(outstem+'_multi.pdf')
plt.close(F)

# Method 2a
#
used_2a = False
if len(tslicesM2a)!=4 or tslicesM2a[-1]<=tslicesM2a[-2]:
  print ('Error: tslicesM2a =',tslicesM2a,'does not have length 4 or has insufficient span.')
  print ('Skipping Method 2a ...')
else:
  # Proceed to implement Method 2a
  used_2a = True
  print ('Alt. time slices (Method 2a): ',tslicesM2a)
  tfmin = tslicesM2a[2]; tfmax = tslicesM2a[3]
  ntM2a = tfmax-tfmin+1
  print ('Method 2a, progress of calculation:')
  sys.stdout.write('|')
  for iy in range(ny): sys.stdout.write(' ')
  print ('| <- 100%')
  sys.stdout.write('|')
  Method2a_slopes = numpy.zeros((ny,nx))
  Method2a_vals = numpy.zeros((ny,nx,ntM2a))
  lngraw = numpy.zeros((ntM2a))
  for iy in range(ny):
    sys.stdout.write('*'); sys.stdout.flush()
    for ix in range(nx):
      if is_good[iy,ix]==1:
        for t in range(ntM2a):
          temp_tslices = [tslicesM2a[0], tslicesM2a[1], tslicesM2a[1], tfmin+t]
          if fullref:
            lightref = pyirc.ref_array_onerow(lightfiles, formatpars, iy, ny, temp_tslices, False)
            darkref = pyirc.ref_array_onerow(darkfiles, formatpars, iy, ny, temp_tslices, False)
          region_cube = pyirc.pixel_data(lightfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], temp_tslices,
                        [sensitivity_spread_cut, True], False)
          dark_cube = pyirc.pixel_data(darkfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], temp_tslices,
                        [sensitivity_spread_cut, False], False)
          info = pyirc.basic(region_cube, dark_cube, temp_tslices, lightref[:,iy,:], darkref[:,iy,:], basicpar, False)
          Method2a_vals[iy,ix,t] = lngraw[t] = numpy.log(info[1])
        # Build least squares fit
        mS, cS = numpy.linalg.lstsq(numpy.vstack([numpy.array(range(ntM2a)), numpy.ones(ntM2a)]).T, lngraw, rcond=-1)[0]
        Method2a_slopes[iy,ix] = mS/full_info[iy,ix,pyirc.swi.I]
  print ('|')
  print ('Mean slope d[ln graw]/d[I td] at fixed ta,tb =', numpy.mean(is_good*Method2a_slopes)/numpy.mean(is_good))
  print ('')
  # Predicted slopes
  slope_2a_BFE = 3*mean_full_info[pyirc.swi.beta] - (1+4*mean_full_info[pyirc.swi.alphaH]
                 +4*mean_full_info[pyirc.swi.alphaV])*mean_full_info[pyirc.swi.ker0]
  slope_2a_NLIPC = 3*mean_full_info[pyirc.swi.beta] - 2*(1+4*mean_full_info[pyirc.swi.alphaH]
                   +4*mean_full_info[pyirc.swi.alphaV])*mean_full_info[pyirc.swi.ker0]

# Method 2b
#
used_2b = False
if len(tslicesM2b)!=4 or tslicesM2b[-1]<=tslicesM2b[-2]:
  print ('Error: tslicesM2b =',tslicesM2b,'does not have length 4 or has insufficient span.')
  print ('Skipping Method 2b ...')
else:
  # Proceed to implement Method 2b
  used_2b = True
  print ('Alt. time slices (Method 2b): ',tslicesM2b)
  tfminB = tslicesM2b[2]; tfmaxB = tslicesM2b[3]
  ntM2b = tfmaxB-tfminB+1
  print ('Method 2b, progress of calculation:')
  sys.stdout.write('|')
  for iy in range(ny): sys.stdout.write(' ')
  print ('| <- 100%')
  sys.stdout.write('|')
  Method2b_slopes = numpy.zeros((ny,nx))
  Method2b_vals = numpy.zeros((ny,nx,ntM2b))
  lngraw = numpy.zeros((ntM2b))
  for iy in range(ny):
    sys.stdout.write('*'); sys.stdout.flush()
    for ix in range(nx):
      if is_good[iy,ix]==1:
        for t in range(ntM2b):
          temp_tslices = [tslicesM2b[0]+t, tslicesM2b[1]+t, tslicesM2b[1]+t, tslicesM2b[2]+t]
          if fullref:
            lightref = pyirc.ref_array_onerow(lightfiles, formatpars, iy, ny, temp_tslices, False)
            darkref = pyirc.ref_array_onerow(darkfiles, formatpars, iy, ny, temp_tslices, False)
          region_cube = pyirc.pixel_data(lightfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], temp_tslices,
                        [sensitivity_spread_cut, True], False)
          dark_cube = pyirc.pixel_data(darkfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], temp_tslices,
                        [sensitivity_spread_cut, False], False)
          info = pyirc.basic(region_cube, dark_cube, temp_tslices, lightref[:,iy,:], darkref[:,iy,:], basicpar, False)
          Method2b_vals[iy,ix,t] = lngraw[t] = numpy.log(info[1])
        # Build least squares fit
        mS, cS = numpy.linalg.lstsq(numpy.vstack([numpy.array(range(ntM2b)), numpy.ones(ntM2b)]).T, lngraw, rcond=-1)[0]
        Method2b_slopes[iy,ix] = mS/full_info[iy,ix,pyirc.swi.I]
  print ('|')
  print ('Mean slope d[ln graw]/d[I tb] at fixed tab,tad =', numpy.mean(is_good*Method2b_slopes)/numpy.mean(is_good))
  print ('')
  # Predicted slopes
  slope_2b_BFE = 2*mean_full_info[pyirc.swi.beta]
  slope_2b_NLIPC = 2*mean_full_info[pyirc.swi.beta] + 2*(1+4*mean_full_info[pyirc.swi.alphaH]
                   +4*mean_full_info[pyirc.swi.alphaV])*mean_full_info[pyirc.swi.ker0]

# Method 3
#
used_3 = False
if len(tslicesM3)!=4 or tslicesM3[-1]<=tslicesM3[-2]:
  print ('Error: tslicesM3 =',tslicesM3,'does not have length 4 or has insufficient span.')
  print ('Skipping Method 3 ...')
else:
  # Proceed to implement Method 3
  used_3 = True
  print ('Alt. time slices (Method 3): ',tslicesM3)
  tfmin3 = tslicesM3[2]; tfmax3 = tslicesM3[3]
  ntM3 = tfmax3-tfmin3+1
  print ('Method 3, progress of calculation:')
  sys.stdout.write('|')
  for iy in range(ny): sys.stdout.write(' ')
  print ('| <- 100%')
  sys.stdout.write('|')
  Method3_slopes = numpy.zeros((ny,nx))
  Method3_vals = numpy.zeros((ny,nx,ntM3))
  Method3_alphas = numpy.zeros((ny,nx,ntM3))
  CCraw = numpy.zeros((ntM3))
  for iy in range(ny):
    sys.stdout.write('*'); sys.stdout.flush()
    for ix in range(nx):
      if is_good[iy,ix]==1:
        for t in range(ntM3):
          temp_tslices = [tslicesM3[0], tslicesM3[1], tslicesM3[0], tfmin3+t]
          if fullref:
            lightref = pyirc.ref_array_onerow(lightfiles, formatpars, iy, ny, temp_tslices, False)
            darkref = pyirc.ref_array_onerow(darkfiles, formatpars, iy, ny, temp_tslices, False)
          region_cube = pyirc.pixel_data(lightfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], temp_tslices,
                        [sensitivity_spread_cut, True], False)
          dark_cube = pyirc.pixel_data(darkfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], temp_tslices,
                        [sensitivity_spread_cut, False], False)
          info = pyirc.basic(region_cube, dark_cube, temp_tslices, lightref[:,iy,:], darkref[:,iy,:], basicpar, False)
          Method3_vals[iy,ix,t] = CCraw[t] = (info[pyirc.swi.tCH]+info[pyirc.swi.tCV])/2.*full_info[iy,ix,pyirc.swi.g]**2\
            /(full_info[iy,ix,pyirc.swi.I]*(temp_tslices[-1]-temp_tslices[0]))
          Method3_alphas[iy,ix,t] = (info[pyirc.swi.alphaH]+info[pyirc.swi.alphaV])/2.
        # Build least squares fit
        mS, cS = numpy.linalg.lstsq(numpy.vstack([numpy.array(range(ntM3)), numpy.ones(ntM3)]).T, CCraw, rcond=-1)[0]
        Method3_slopes[iy,ix] = mS/full_info[iy,ix,pyirc.swi.I]
  print ('|')
  print ('Mean slope d[g^2/(Itad) Cadj,ad]/d[I td] at fixed ta,tb =', numpy.mean(is_good*Method3_slopes)/numpy.mean(is_good))
  print ('')
  # Predicted slopes
  ave = (mean_full_info[pyirc.swi.ker0-1]+mean_full_info[pyirc.swi.ker0+1]+mean_full_info[pyirc.swi.ker0-(2*pyirc.swi.s+1)]
        +mean_full_info[pyirc.swi.ker0+(2*pyirc.swi.s+1)])/4.
  slope_3_beta = -4*(mean_full_info[pyirc.swi.alphaH]+mean_full_info[pyirc.swi.alphaV])*mean_full_info[pyirc.swi.beta]
  slope_3_BFE = -4*(mean_full_info[pyirc.swi.alphaH]+mean_full_info[pyirc.swi.alphaV])*mean_full_info[pyirc.swi.beta] + ave
  slope_3_NLIPC = -4*(mean_full_info[pyirc.swi.alphaH]+mean_full_info[pyirc.swi.alphaV])*mean_full_info[pyirc.swi.beta] + ave*2.

# Non-linearity corrections, Methods 2 and 3:
#
# I have commented out the outputs for this section. We may put them back as needed.
#for it in range(ntSub): 
#  print ('{:2d} {:8.2f} {:8.2f} {:8.2f}'.format(it+1, nlMean[it], nlFit[it], nlMeanDer[it]))
print ('Non-linearity correction tables:')
if used_2a:
  print ('2a:')
  vec = []
  for t in range(tslicesM2a[2], tslicesM2a[3]+1):
    offsets = pyirc.compute_gain_corr_many(nlfit, nlder, full_info[:,:,pyirc.swi.I]*full_info[:,:,pyirc.swi.beta], [tslicesM2a[0],tslicesM2a[1],t], basicpar.reset_frame, is_good)
    print (t, numpy.mean(offsets*is_good)/numpy.mean(is_good))
    vec += [numpy.mean(offsets*is_good)/numpy.mean(is_good)]
  PV2a = max(vec)-min(vec); print ('PV: ', PV2a)
if used_2b:
  print ('2b:')
  vec = []
  dt1 = tslicesM2b[1] - tslicesM2b[0]
  dt2 = tslicesM2b[2] - tslicesM2b[0]
  for t in range(tslicesM2b[0], tslicesM2b[3]-tslicesM2b[2]+1):
    offsets = pyirc.compute_gain_corr_many(nlfit, nlder, full_info[:,:,pyirc.swi.I]*full_info[:,:,pyirc.swi.beta], [t,t+dt1,t+dt2], basicpar.reset_frame, is_good)
    print (t, numpy.mean(offsets*is_good)/numpy.mean(is_good))
    vec += [numpy.mean(offsets*is_good)/numpy.mean(is_good)]
  PV2b = max(vec)-min(vec); print ('PV: ', PV2b)
if used_3:
  print ('3:')
  vec = []
  for t in range(tslicesM3[2], tslicesM3[3]+1):
    offsets = pyirc.compute_xc_corr_many(nlfit, nlder, full_info[:,:,pyirc.swi.I]*full_info[:,:,pyirc.swi.beta], [tslicesM3[0],t], basicpar.reset_frame, is_good)
    alpha3 = (full_info[:,:,pyirc.swi.alphaH]+full_info[:,:,pyirc.swi.alphaV])/2.
    offsets *= 2. * alpha3 * (1.-4*alpha3)
    print (t, numpy.mean(offsets*is_good)/numpy.mean(is_good))
    vec += [numpy.mean(offsets*is_good)/numpy.mean(is_good)]
  PV3 = max(vec)-min(vec); print ('PV: ', PV3)
print ('')

# Method 2 and 3 characterization
# Multi-panel figure showing basic characterization
matplotlib.rcParams.update({'font.size': 8})
F = plt.figure(figsize=(3.5,9))
if used_2a:
  S = F.add_subplot(3,1,1)
  S.set_title(r'Raw gain vs. interval duration')
  S.set_xlabel(r'Signal level $It_{'+'{:d}'.format(tslicesM2a[0])+r',d}$ [ke]')
  S.set_ylabel(r'$\ln g^{\rm raw}_{' +'{:d},{:d}'.format(tslicesM2a[0],tslicesM2a[1]) +r',d}$')
  SX = [numpy.mean(is_good*full_info[:,:,pyirc.swi.I]*myt)/numpy.mean(is_good)/1.0e3 for myt in range(tfmin-tslicesM2a[0], tfmax+1-tslicesM2a[0])]
  SY = [numpy.mean(is_good*Method2a_vals[:,:,t])/numpy.mean(is_good) for t in range(ntM2a)]
  SS = [] # std. dev. on the mean
  for t in range(ntM2a):
    SS += [ numpy.sqrt((numpy.mean(is_good*Method2a_vals[:,:,t]**2)/numpy.mean(is_good)-SY[t]**2)/(numpy.sum(is_good)-1)) ]
  xc = numpy.mean(numpy.array(SX))
  yc = numpy.mean(numpy.array(SY))
  S.set_xlim(min(SX)-.05*(max(SX)-min(SX)), max(SX)+.05*(max(SX)-min(SX)))
  xr = numpy.arange(min(SX), max(SX), (max(SX)-min(SX))/256.)
  S.errorbar([xc], [min(SY)], yerr=[PV2a/2.], marker=',', color='k', ls='None')
  S.text(xc+.05*(max(SX)-min(SX)), min(SY), 'sys nl', color='k')
  S.errorbar(SX, SY, yerr=SS, marker='x', color='r', ls='None')
  S.plot(xr, yc+(xr-xc)*slope_2a_BFE*1e3, 'g--', label='pure BFE')
  S.plot(xr, yc+(xr-xc)*slope_2a_NLIPC*1e3, 'b-', label='pure NL-IPC')
  S.legend(loc=2)
if used_2b:
  S = F.add_subplot(3,1,2)
  S.set_title(r'Raw gain vs. interval center')
  S.set_xlabel(r'Signal level $It_{'+'{:d}'.format(tslicesM2b[0])+r',a}$ [ke]')
  S.set_ylabel(r'$\ln g^{\rm raw}_{' +'a,a+{:d},a+{:d}'.format(tslicesM2b[1]-tslicesM2b[0],tslicesM2b[2]-tslicesM2b[0]) +r'}$')
  SX = [numpy.mean(is_good*full_info[:,:,pyirc.swi.I]*myt)/numpy.mean(is_good)/1.0e3 for myt in range(ntM2b)]
    # the -1e-5 is to set the x-axis and has no effect
  SY = [numpy.mean(is_good*Method2b_vals[:,:,t])/numpy.mean(is_good) for t in range(ntM2b)]
  SS = [] # std. dev. on the mean
  for t in range(ntM2b):
    SS += [ numpy.sqrt((numpy.mean(is_good*Method2b_vals[:,:,t]**2)/numpy.mean(is_good)-SY[t]**2)/(numpy.sum(is_good)-1)) ]
  xc = numpy.mean(numpy.array(SX))
  yc = numpy.mean(numpy.array(SY))
  S.set_xlim(min(SX)-.05*(max(SX)-min(SX)), max(SX)+.05*(max(SX)-min(SX)))
  xr = numpy.arange(min(SX), max(SX), (max(SX)-min(SX))/256.)
  S.errorbar([xc], [min(SY)], yerr=[PV2b/2.], marker=',', color='k', ls='None')
  S.text(xc+.05*(max(SX)-min(SX)), min(SY), 'sys nl', color='k')
  S.errorbar(SX, SY, yerr=SS, marker='x', color='r', ls='None')
  S.plot(xr, yc+(xr-xc)*slope_2b_BFE*1e3, 'g--', label='pure BFE')
  S.plot(xr, yc+(xr-xc)*slope_2b_NLIPC*1e3, 'b-', label='pure NL-IPC')
  S.legend(loc=2)
if used_3:
  S = F.add_subplot(3,1,3)
  S.set_title(r'CDS ACF vs. signal')
  S.set_xlabel(r'Signal level $It_{'+'{:d}'.format(tslicesM3[0])+r',d}$ [ke]')
  S.set_ylabel(r'$g^2C_{'+'{:d}'.format(tslicesM3[0])+r'd'+'{:d}'.format(tslicesM3[0])+r'd}(\langle1,0\rangle)/[It_{'\
    +'{:d}'.format(tslicesM3[0])+r'd}]$')
  SX = [numpy.mean(is_good*full_info[:,:,pyirc.swi.I]*myt)/numpy.mean(is_good)/1.0e3 for myt in range(tfmin3-tslicesM3[0], tfmax3+1-tslicesM3[0])]
  SY = [numpy.mean(is_good*Method3_vals[:,:,t])/numpy.mean(is_good) for t in range(ntM3)]
  SS = [] # std. dev. on the mean
  for t in range(ntM3):
    SS += [ numpy.sqrt((numpy.mean(is_good*Method3_vals[:,:,t]**2)/numpy.mean(is_good)-SY[t]**2)/(numpy.sum(is_good)-1)) ]
  xc = numpy.mean(numpy.array(SX))
  yc = numpy.mean(numpy.array(SY))
  S.set_xlim(min(SX)-.05*(max(SX)-min(SX)), max(SX)+.05*(max(SX)-min(SX)))
  xr = numpy.arange(min(SX), max(SX), (max(SX)-min(SX))/256.)
  S.errorbar([xc], [min(SY)], yerr=[PV3], marker=',', color='k', ls='None')
  S.text(xc+.05*(max(SX)-min(SX)), min(SY)+PV3, 'sys nl', color='k')
  S.errorbar(SX, SY, yerr=SS, marker='x', color='r', ls='None')
  S.plot(xr, yc+(xr-xc)*slope_3_BFE*1e3, 'g--', label='pure BFE')
  S.plot(xr, yc+(xr-xc)*slope_3_NLIPC*1e3, 'b-', label='pure NL-IPC')
  S.plot(xr, yc+(xr-xc)*slope_3_beta*1e3, 'k:', label='beta only')
  S.legend(loc=2)
  print ('Method 3 implied slopes =', slope_3_beta, slope_3_BFE, slope_3_NLIPC)
F.set_tight_layout(True)
F.savefig(outstem+'_m23.pdf')
plt.close(F)

# Text output
thisOut = open(outstem+'_summary.txt', 'w')
# Print information in the file header
thisOut.write('# This summary created at {:s}\n'.format(time.asctime(time.localtime(time.time()))))
thisOut.write('# Uses pyirc v{:d}\n'.format(pyirc.get_version()))
thisOut.write('# Detector: '+mydet+'\n')
thisOut.write('#\n# Files used:\n')
thisOut.write('# Light:\n')
for f in lightfiles: thisOut.write('#   {:s}\n'.format(f))
thisOut.write('# Dark:\n')
for f in darkfiles: thisOut.write('#   {:s}\n'.format(f))
thisOut.write('# Input format {:d}\n'.format(formatpars))
thisOut.write('# Time slices:')
for t in tslices: thisOut.write(' {:3d}'.format(t))
thisOut.write('\n')
thisOut.write('# Mask: ' + str(maskX) + ',' + str(maskY) + '\n')
thisOut.write('#\n')
thisOut.write('# Cut on good pixels {:7.4f}% deviation from median\n'.format(100*sensitivity_spread_cut))
thisOut.write('# Dimensions: {:3d}(x) x {:3d}(y) super-pixels, {:4d} good\n'.format(nx,ny,int(numpy.sum(is_good))))
thisOut.write('# Frame number corresponding to reset: {:d}\n'.format(basicpar.reset_frame))
thisOut.write('# Reference pixel subtraction for linearity: {:s}\n'.format(str(fullref)))
thisOut.write('# Quantile for variance computation = {:9.6f}%\n'.format(basicpar.g_ptile))
thisOut.write('# Clipping fraction epsilon = {:9.7f}\n'.format(basicpar.epsilon))
thisOut.write('# Lead-trail subtraction for IPC correlation = ' + str(basicpar.leadtrailSub) + '\n')
thisOut.write('# Characterization type: '+mychar+'\n')
if mychar.lower()=='advanced':
  thisOut.write('#   dt = {:d},{:d}, ncycle = {:d}\n'.format(tchar1,tchar2,ncycle))
thisOut.write('# Non-linearity settings: basicpar.fullnl={:s} bfepar.fullnl={:s} basicpar.use_allorder={:s}\n'.format(
  str(basicpar.fullnl), str(bfepar.fullnl), str(basicpar.use_allorder) ))
thisOut.write('# BFE Method 1\n#   Baseline subtraction = {:s}\n'.format(str(bfepar.blsub)))
thisOut.write('# BFE Method 2a\n#   Enabled = {:s}\n'.format(str(used_2a)))
thisOut.write('# BFE Method 2b\n#   Enabled = {:s}\n'.format(str(used_2b)))
thisOut.write('# BFE Method 3\n#   Enabled = {:s}\n'.format(str(used_3)))
thisOut.write('# Hot pixels\n#   Enabled = {:s}\n'.format(str(hotpix)))
if hotpix:
  thisOut.write('#   Parameters = {:s}\n'.format(str(hotpix_ADU_range)))
  if ref_for_hotpix_is_autocorr:
    thisOut.write('#   ref for delta alpha = autocorr\n')
  else:
    thisOut.write('#   ref for delta alpha = last frame used\n')
  if hotpix_slidemed:
    thisOut.write('#   median method = sliding\n')
  else:
    thisOut.write('#   median method = normal\n')
thisOut.write('# Associated figures:\n')
thisOut.write('#   {:s}\n'.format(outstem+'_multi.pdf'))
thisOut.write('#   {:s}\n'.format(outstem+'_m23.pdf'))
thisOut.write('#   {:s}\n'.format(outstem+'_hotipc.pdf'))
thisOut.write('#\n')
thisOut.write('# Columns:\n'); col=1
thisOut.write('# {:3d}, X (super pixel grid)\n'.format(col)); col+=1
thisOut.write('# {:3d}, Y (super pixel grid)\n'.format(col)); col+=1
thisOut.write('# {:3d}, number of good pixels\n'.format(col)); col+=1
thisOut.write('# {:3d}, raw gain (e/DN)\n'.format(col)); col+=1
thisOut.write('# {:3d}, alpha-corrected gain (e/DN)\n'.format(col)); col+=1
thisOut.write('# {:3d}, alpha,beta-corrected gain (e/DN)\n'.format(col)); col+=1
thisOut.write('# {:3d}, IPC alpha horizontal\n'.format(col)); col+=1
thisOut.write('# {:3d}, IPC alpha vertical\n'.format(col)); col+=1
thisOut.write('# {:3d}, nonlinearity beta (e^-1)\n'.format(col)); col+=1
thisOut.write('# {:3d}, charge per time slice (e)\n'.format(col)); col+=1
thisOut.write('# {:3d}, IPC alpha diagonal (if computed; otherwise all 0s)\n'.format(col)); col+=1
thisOut.write('# {:3d}, C_H at slices {:d},{:d} (DN^2)\n'.format(col, tslices[0], tslices[-1])); col+=1
thisOut.write('# {:3d}, C_V at slices {:d},{:d} (DN^2)\n'.format(col, tslices[0], tslices[-1])); col+=1
for jb in range(2*pyirc.swi.s+1):
  for ib in range(2*pyirc.swi.s+1):
    thisOut.write('# {:3d}, BFE kernel K^2a (+NL-IPC) at ({:2d},{:2d}) (e^-1)\n'.format(col, ib-pyirc.swi.s, jb-pyirc.swi.s)); col+=1
if pyirc.swi.p>0:
  thisOut.write('# {:3d}, time intercept\n'.format(col)); col += 1
  for ip in range(2, pyirc.swi.p+1):
    thisOut.write('# {:3d}, additional non-linearity coefficient, order {:d} (DN^-{:d})\n'.format(col, ip, ip-1)); col+=1
if used_2a: thisOut.write('# {:3d}, Method 2a slope (e^-1)\n'.format(col)); col+=1
if used_2b: thisOut.write('# {:3d}, Method 2b slope (e^-1)\n'.format(col)); col+=1
if used_3: thisOut.write('# {:3d}, Method 3 slope (e^-1)\n'.format(col)); col+=1
thisOut.write('#\n')
# Now make the data table
for iy in range(ny):
  for ix in range(nx):
    # Print the column first, then row (normal human-read order, note this is the reverse of internal Python)
    thisOut.write('{:3d} {:3d}'.format(ix,iy))
    for col in range(my_dim): thisOut.write(' {:14.7E}'.format(full_info[iy,ix,col]))
    if used_2a: thisOut.write(' {:14.7E}'.format(Method2a_slopes[iy,ix]))
    if used_2b: thisOut.write(' {:14.7E}'.format(Method2b_slopes[iy,ix]))
    if used_3: thisOut.write(' {:14.7E}'.format(Method3_slopes[iy,ix]))
    thisOut.write('\n')
thisOut.close()

if hotpix:
  print ('Start hot pixels ...')
  hotY, hotX = pyirc.hotpix(darkfiles, formatpars, range(1,NTMAX), hotpix_ADU_range, True)
  print ('Number of pixels selected:', len(hotX)) # only printed for de-bugging -> , len(hotY)
  dtstep = 5 # <-- right now this is hard coded
  htsteps = range(1,NTMAX,dtstep)
  if hotpix_logtspace:
    htsteps = [1]
    for k in range(1,12):
      if 2**k<NTMAX-1: htsteps += [2**k]
      if k>=2:
        if 5*2**(k-2)<NTMAX-1: htsteps += [5*2**(k-2)]
      if 3*2**(k-1)<NTMAX-1: htsteps += [3*2**(k-1)]
      if k>=2:
        if 7*2**(k-2)<NTMAX-1: htsteps += [7*2**(k-2)]
    htsteps += [NTMAX-1]
  beta_gain = full_info[:,:,pyirc.swi.beta]*full_info[:,:,pyirc.swi.g]
  print (beta_gain)
  hotcube = pyirc.hotpix_ipc(hotY, hotX, darkfiles, formatpars, htsteps, [beta_gain, False], True)
  nhstep = len(htsteps)
  print ('number of time steps ->', nhstep)
  fromcorr_alpha = numpy.zeros((len(hotX)))
  hotpix_alpha = numpy.zeros((len(hotX), nhstep))
  hotpix_alpha_num = numpy.zeros((len(hotX), nhstep))
  hotpix_alpha_den = numpy.zeros((len(hotX), nhstep))
  hotpix_alphaD = numpy.zeros((len(hotX), nhstep))
  hotpix_signal = numpy.zeros((len(hotX), nhstep))
  #
  # generate and write hot pixel information
  thisOut = open(outstem+'_hot.txt', 'w')
  for jpix in range(len(hotX)):
    iy = hotY[jpix]//dy
    ix = hotX[jpix]//dx
    fromcorr_alpha[jpix] = full_info[iy,ix,pyirc.swi.alphaH]/2.+full_info[iy,ix,pyirc.swi.alphaV]/2.
    thisOut.write('{:4d} {:4d} {:8.6f}'.format(hotX[jpix], hotY[jpix], fromcorr_alpha[jpix]))
    for t in range(1,nhstep):
      R = ( numpy.mean(hotcube[jpix,t,1:5]) - hotcube[jpix,t,-1] ) / (hotcube[jpix,t,0]-hotcube[jpix,t,-1] )
      S = ( numpy.mean(hotcube[jpix,t,5:9]) - hotcube[jpix,t,-1] ) / (hotcube[jpix,t,0]-hotcube[jpix,t,-1] )
      hotpix_alpha[jpix, t] = R/(1.+4*R+4*S)
      hotpix_alpha_num[jpix, t] = R
      hotpix_alpha_den[jpix, t] = 1.+4*R+4*S
      hotpix_alphaD[jpix, t] = S/(1.+4*R+4*S)
      hotpix_signal[jpix, t] = hotcube[jpix,t,0]-hotcube[jpix,t,-1]
      thisOut.write(' {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}'.format(hotcube[jpix,t,0], numpy.mean(hotcube[jpix,t,1:5]), hotcube[jpix,t,-1],
        (hotcube[jpix,t,1]+hotcube[jpix,t,3]-hotcube[jpix,t,2]-hotcube[jpix,t,4])/4., numpy.mean(hotcube[jpix,t,5:9])))
    thisOut.write('\n')
  thisOut.close()

  # report median levels
  print ('IPC relative to nominal (signal, median, uncert):')
  ipcmed_x = numpy.zeros((nhstep))
  ipcmed_y = numpy.zeros((nhstep))
  ipcmed_yerr = numpy.zeros((nhstep))
  delta = .5/numpy.sqrt(len(hotX))
  for t in range(1,nhstep):
    my_y = hotpix_alpha[:,t]-hotpix_alpha[:,-1]
    if ref_for_hotpix_is_autocorr: my_y = hotpix_alpha[:,t] - fromcorr_alpha
    ipcmed_x[t] = numpy.nanpercentile(hotpix_signal[:, t], 50.)
    if hotpix_slidemed:
      X_ = hotpix_alpha_den[:,t]
      Y_ = hotpix_alpha_num[:,t]-fromcorr_alpha*hotpix_alpha_den[:,t]
      ipcmed_y[t] = pyirc.slidemed_percentile(X_, Y_, 50)
      ipcmed_yerr[t] = ( pyirc.slidemed_percentile(X_, Y_, 50.+100*delta)
        - pyirc.slidemed_percentile(X_, Y_, 50.-100*delta) ) / 2.
    else:
      ipcmed_y[t] = numpy.nanpercentile(my_y, 50.)
      ipcmed_yerr[t] = (numpy.nanpercentile(my_y, 50.+100*delta)-numpy.nanpercentile(my_y, 50.-100*delta))/2.
    print ('{:10.2f} {:9.6f} {:9.6f}'.format(ipcmed_x[t], ipcmed_y[t], ipcmed_yerr[t]))
  print ('')
  print ('median alphaD:', '{:9.6f} {:9.6f}'.format(numpy.nanpercentile(hotpix_alphaD[:,-1], 50.),
    (numpy.nanpercentile(hotpix_alphaD[:,-1], 50.+100*delta)-numpy.nanpercentile(hotpix_alphaD[:,-1], 50.-100*delta))/2.))

  # bigger grid for IPC comparisons
  NG=4
  grid_alphaCorr = numpy.zeros((NG,NG)); grid_alphaCorrErr = numpy.zeros((NG,NG))
  grid_alphaHot = numpy.zeros((NG,NG)); grid_alphaHotErr = numpy.zeros((NG,NG))
  if ny%NG==0 and nx%NG==0:
    stepx = nx//NG; stepy = ny//NG
    sp = N//NG;
    for ix in range(NG):
      for iy in range(NG):
        # bin the auto-correlation measurements
        suba = full_info[iy*stepy:(iy+1)*stepy, ix*stepx:(ix+1)*stepx, :]
        mya = (suba[:,:,pyirc.swi.alphaH]+suba[:,:,pyirc.swi.alphaV])/2.
        pmask = suba[:,:,0] > 0
        n = mya[pmask].size
        if n>1:
          grid_alphaCorr[iy,ix] = mya[pmask].mean()
          grid_alphaCorrErr[iy,ix] = mya[pmask].std()/numpy.sqrt(n-1)
        #
        # bin the hot pixel measurements -- use final time slice!
        u = hotpix_alpha[numpy.logical_and(hotY//sp==iy, hotX//sp==ix), -1]
        if u.size>1:
          grid_alphaHot[iy,ix] = numpy.nanpercentile(u,50)
          delta = .5/numpy.sqrt(u.size)
          grid_alphaHotErr[iy,ix] = (numpy.nanpercentile(u, 50.+100*delta)-numpy.nanpercentile(u, 50.-100*delta))/2.
  print (grid_alphaCorr, grid_alphaCorrErr, grid_alphaHot, grid_alphaHotErr)

  # hot pixel plots
  matplotlib.rcParams.update({'font.size': 8})
  if narrowfig:
    F = plt.figure(figsize=(3.5,6))
  else:
    F = plt.figure(figsize=(7,6))
  #
  # hot pixel locations
  if narrowfig:
    S = F.add_subplot(2,1,1)
  else:
    S = F.add_subplot(2,2,1)
  S.set_title('hot pixel locations: '+mydet)
  S.set_xlim(0,N); S.set_ylim(0,N); S.set_aspect('equal')
  S.xaxis.set_ticks(numpy.linspace(0,N,num=5)); S.yaxis.set_ticks(numpy.linspace(0,N,num=5))
  S.grid(True, color='g', linestyle='-')
  S.scatter(hotX+.5, hotY+.5, s=3, marker='.', color='r')
  #
  # these two panels only in the full version, not the narrow version
  if not narrowfig:
    # hot pixel level
    S = F.add_subplot(2,2,2)
    S.set_xlabel(r'Signal level $S_{1,' + '{:d}'.format(htsteps[-1]) +'}$ [DN]')
    S.set_ylabel(r'IPC $\alpha$ [%]')
    SX = hotpix_signal[:,-1]
    SY = hotpix_alpha[:,-1]/.01
    S.set_title(r'IPC $\alpha$ for hot pixels')
    S.set_xlim(.95*(htsteps[-1]-1)/(NTMAX-1.0)*hotpix_ADU_range[0], 1.05*hotpix_ADU_range[1])
    S.set_ylim(0,4.)
    S.grid(True, color='g', linestyle='-')
    S.scatter(SX, SY, s=3, marker='.', color='r')
    #
    # dependence on signal level
    S = F.add_subplot(2,2,3)
    S.set_xlabel(r'Signal level $S_{1,b}$ [DN]')
    S.set_ylabel(r'IPC $\alpha(S_{1,b})-\alpha(S_{1,' + '{:d}'.format(htsteps[-1]) + '})$ [%]')
    for t in range(1,nhstep):
      SXa = hotpix_signal[:,t]
      SYa = (hotpix_alpha[:,t]-hotpix_alpha[:,-1])/.01
      if t==1:
        SX = SXa; SY = SYa
      else:
        print (t, numpy.shape(SX), numpy.shape(SXa))
        SX = numpy.concatenate((SX,SXa)); SY = numpy.concatenate((SY,SYa))
    S.set_title(r'IPC signal dependence $\Delta\alpha$')
    S.set_xlim(0., hotpix_ADU_range[1])
    S.set_ylim(-1.5,1.5)
    S.xaxis.set_ticks(numpy.linspace(0,hotpix_ADU_range[1],num=6)); S.yaxis.set_ticks(numpy.linspace(-1.5,1.5,num=11))
    S.grid(True, color='g', linestyle='-', linewidth=.5)
    S.scatter(SX, SY, s=.25, marker='+', color='r')
    S.errorbar(ipcmed_x[1:], ipcmed_y[1:]/.01, yerr=ipcmed_yerr[1:]/.01, ms=2, marker='o', color='k', ls='None')
  #
  # comparison with auto-correlations
  if narrowfig:
    S = F.add_subplot(2,1,2)
  else:
    S = F.add_subplot(2,2,4)
  S.set_title(r'hot pixels vs. autocorr. IPC $\alpha$')
  scale_test = numpy.concatenate((grid_alphaCorr.flatten(), grid_alphaHot.flatten()))
  smin = 0.92 * numpy.min(scale_test[scale_test>0]) / .01
  smax = 1.08 * numpy.max(scale_test) / .01
  S.set_xlim(smin,smax); S.set_ylim(smin,smax); S.set_aspect('equal')
  S.set_xlabel(r'autocorrelation $\alpha$ [%]'); S.set_ylabel(r'hot pixel $\alpha$ [%]')
  S.grid(True, color='g', linestyle='-', linewidth=.5)
  S.errorbar(grid_alphaCorr.flatten()/.01, grid_alphaHot.flatten()/.01,
    xerr=grid_alphaCorrErr.flatten()/.01, yerr=grid_alphaHotErr.flatten()/.01,
    ms=1, marker='o', color='k', capsize=1, ls='None')
  xr = numpy.linspace(0,4,num=65)
  S.plot(xr, xr, 'r-')

  F.set_tight_layout(True)
  F.savefig(outstem+'_hotipc.pdf')
  plt.close(F)

