import os
import sys
import time
import re
import numpy
import pyirc
import ftsolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import copy

class EmptyClass:
  pass

outstem = 'default_output'
use_cmap = 'gnuplot'

mydet = ''
lightfiles = []
darkfiles = []
vislightfiles = []
visdarkfiles = []
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
basicpar.vis_med_correct = False

# Parameters for BFE
bfepar = EmptyClass()
bfepar.epsilon = .01
bfepar.treset = basicpar.reset_frame
bfepar.blsub = True
bfepar.fullnl = False
bfepar.vis = True
#
copy_ir_bfe = False

# Plotting parameters
narrowfig = False

# Separate parameters for visible BFE?
has_visbfe = False

# Read in information
config_file = sys.argv[1]
if len(sys.argv)>2:
  ir_flag = sys.argv[2]  # Currently not doing anything with this but could
  print('Running test_run.py without charge diffusion')
  cmd='python test_run.py %s'%config_file
  os.system(cmd)
  #sys.exit()  # Probably we wil always want to continue.
  # If this is run, we may want to read in resulting summary.txt files
  # This is simply duplicating some commands below
  #with open(conf, 'r') as ifile: conf_content=conf.read().splitlines()
  #for line in content:
  #  m=re.search(r'^[A-Z]+\:', line)
  #  m = re.search(r'^OUTPUT\:\s*(\S*)', line)
  #  if m: outstem = m.group(1)
  # This part still being written
    

with open(config_file) as myf: content = myf.read().splitlines()
is_in_light = is_in_dark = is_in_vislight = is_in_visdark = False
maskX = [] # list of regions to mask
maskY = []
for line in content:
  # Cancellations
  m = re.search(r'^[A-Z]+\:', line)
  if m: is_in_light = is_in_dark = is_in_vislight = is_in_visdark = False

  # Searches for files -- must be first given the structure of this script!
  # The visible flats and darks must come after IR flats and darks
  if is_in_light:
    m = re.search(r'^\s*(\S.*)$', line)
    if m: lightfiles += [m.group(1)]
  if is_in_dark:
    m = re.search(r'^\s*(\S.*)$', line)
    if m: darkfiles += [m.group(1)]
  if is_in_vislight:
    m = re.search(r'^\s*(\S.*)$', line)
    if m: vislightfiles += [m.group(1)]
  if is_in_visdark:
    m = re.search(r'^\s*(\S.*)$', line)
    if m: visdarkfiles += [m.group(1)]
        
  # -- Keywords go below here --

  # Search for outputs
  m = re.search(r'^OUTPUT\:\s*(\S*)', line)
  if m: outstem = m.group(1)
  # Search for input files
  m = re.search(r'^LIGHT\:', line)
  if m: is_in_light = True
  m = re.search(r'^DARK\:', line)
  if m: is_in_dark = True
  m = re.search(r'^VISLIGHT\:', line)
  if m: is_in_vislight = True
  m = re.search(r'^VISDARK\:', line)
  if m: is_in_visdark = True

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

  # Visible time stamp range
  m = re.search(r'^VISTIME:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
  if m:
     ts_vis = int(m.group(1))
     te_vis = int(m.group(2))
     tchar1_vis = int(m.group(3))
     tchar2_vis = int(m.group(4))
  #
  m = re.search(r'^VISMEDCORR', line)
  if m: basicpar.vis_med_correct = True
  #
  # Visible BFE
  m = re.search(r'^VISBFETIME:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
  if m:
    tslices_visbfe = [ int(m.group(x)) for x in range(1,5)]
    has_visbfe = True

  m = re.search(r'^COPYIRBFE', line)
  if m: copy_ir_bfe = True

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

# copy to output file
os.system('cp ' + config_file + ' ' + outstem + '_config.txt')

# replace visible time slices for BFE
if has_visbfe: tslices = tslices_visbfe

# set up array size parameters
pyirc.swi.addbfe(s_bfe)
pyirc.swi.addhnl(p_order)
print ('Number of output field per superpixel =', pyirc.swi.N)

# Check number of slices available
NTMAX = 16384
for f in lightfiles+darkfiles:
  nt = pyirc.get_num_slices(formatpars, f)
  if nt<NTMAX: NTMAX=nt

# Copy basicpar parameters to bfepar
bfepar.use_allorder = basicpar.use_allorder

print ('Output will be directed to {:s}*'.format(outstem))
print ('Light files:', lightfiles)
print ('Dark files:', darkfiles)
print ('Visible light files:', vislightfiles)
print ('"Visible" dark files:', visdarkfiles)
print ('Time slices:', tslices, 'max=',NTMAX)
print ('Mask regions:', maskX, maskY)
# 
if len(lightfiles)!=len(darkfiles) or len(lightfiles)<2:
  print ('Failed: {:d} light files and {:d} dark files'.format(len(lightfiles), len(darkfiles)))
  exit()
if len(vislightfiles)!=len(visdarkfiles) or len(vislightfiles)<2:
  print ('Failed: {:d} visible light files and {:d} visible dark files'.format(len(vislightfiles), len(visdarkfiles)))
  exit()

# Additional parameters
# Size of a block
N = pyirc.get_nside(formatpars)
# Side lengths
dx = N//nx
dy = N//ny
# Pixels in a block
npix = dx*dy

# reference pixel subtraction flag
basicpar.subtr_href = fullref

# more allocations
my_dim = pyirc.swi.N
full_info = numpy.zeros((ny,nx,my_dim))
is_good = numpy.zeros((ny,nx))

info_from_ir = numpy.loadtxt(outstem+'_summary.txt')
for j in range(my_dim): full_info[:,:,j] = info_from_ir[:,j+2].reshape((ny,nx))
is_good = numpy.where(full_info[:,:,pyirc.swi.g]>1e-49, 1, 0)

print('Number of good regions =', numpy.sum(is_good))
print('Lower-left corner ->', full_info[0,0,:])

if p_order==0:
  print('Error: did not include polynomial order')
  exit()

# Get Ie
Ie = numpy.zeros((ny,nx))
Ie_alt = numpy.zeros((ny,nx))
Ie_alt2 = numpy.zeros((ny,nx))

print('computing Ie using', ts_vis, te_vis)
nlcubeX, nlfitX, nlderX, pcoefX = pyirc.gen_nl_cube(
  vislightfiles, formatpars, [basicpar.reset_frame, ts_vis, te_vis], [ny,nx],
  full_info[:,:,0], 'abs', False)
for iy in range(ny):
  for ix in range(nx):
    if pcoefX[1,iy,ix]!=0:
      t = numpy.linspace(ts_vis-basicpar.reset_frame, te_vis-basicpar.reset_frame, te_vis-ts_vis+1)
      Signal = numpy.zeros((te_vis-ts_vis+1))
      for ae in range(pyirc.swi.p+1): Signal += pcoefX[ae,iy,ix]*t**ae
      # iterative NL correction
      LinSignal = numpy.copy(Signal)
      for k in range(32):
        LS2 = numpy.copy(LinSignal)
        LinSignal = numpy.copy(Signal)
        LS2 += (LinSignal[-1]-LinSignal[0])/(te_vis-ts_vis) * (ts_vis-basicpar.reset_frame)
        for o in range(2,pyirc.swi.p+1): LinSignal -= full_info[iy,ix,pyirc.swi.Nbb+o-1]*LS2**o
      Ie[iy,ix] = pcoefX[1,iy,ix] * full_info[iy,ix,pyirc.swi.g]
      Ie_alt[iy,ix] = (LinSignal[-1]-LinSignal[0])/(te_vis-ts_vis) * full_info[iy,ix,pyirc.swi.g]
      Sab = Signal[-1]-Signal[0]
      Ie_alt2[iy,ix] = full_info[iy,ix,pyirc.swi.g]*Sab/(te_vis-ts_vis)
      beta_in_e = -full_info[iy,ix,pyirc.swi.Nbb+1:pyirc.swi.Nbb+pyirc.swi.p]/full_info[iy,ix,pyirc.swi.g]**numpy.linspace(1,pyirc.swi.p-1,num=pyirc.swi.p-1) # in e , -
      for k in range(32):
        btcorr = 0
        for j in range(2,pyirc.swi.p+1): btcorr += beta_in_e[j-2]*Ie_alt2[iy,ix]**(j-1)*(t[-1]**j-t[0]**j)
        Ie_alt2[iy,ix] = full_info[iy,ix,pyirc.swi.g]*Sab/(te_vis-ts_vis-btcorr)
    else:
      is_good[iy,ix] = 0 # error

# we use the alt2 method
Ie[:,:] = Ie_alt2

# get vis:IR Ie ratio information
vis_ir_ratio = Ie/full_info[:,:,pyirc.swi.I]
vis_ir_ratio_good = vis_ir_ratio[is_good>.5]
print('VIS:IR ratio information: ', numpy.shape(vis_ir_ratio_good))
print('min, max =', numpy.amin(vis_ir_ratio_good), numpy.amax(vis_ir_ratio_good))
print('percentiles (5th,50th,95th)', numpy.percentile(vis_ir_ratio_good, 5), numpy.percentile(vis_ir_ratio_good, 50),
  numpy.percentile(vis_ir_ratio_good, 95))
print('')

# Allocate space for visible information
vis_bfek = numpy.zeros((ny,nx,5,5))
vis_Phi = numpy.zeros((ny,nx,5,5))
# omega and charge diffusion covariance
QYomega = numpy.zeros((ny,nx))
cdCov = numpy.zeros((ny,nx,3))
cdNiter = numpy.zeros((ny,nx))

# Get correlation functions in each block
nvis = te_vis - ts_vis - tchar2_vis + 1
print ('Visible flat correlation functions, progress of calculation:')
sys.stdout.write('|')
for iy in range(ny): sys.stdout.write(' ')
print ('| <- 100%')
sys.stdout.write('|')
for iy in range(ny):
  sys.stdout.write('*'); sys.stdout.flush()
  if fullref:
    tslices0 = numpy.asarray([ts_vis, ts_vis+tchar1_vis, ts_vis+tchar2_vis])
    lightref_array = []
    darkref_array = []
    for k in range(nvis):
     tslicesk = (tslices0+k).tolist()
     lightref_array.append(pyirc.ref_array(vislightfiles, formatpars, ny, tslicesk, False))
     darkref_array.append(pyirc.ref_array(vislightfiles, formatpars, ny, tslicesk, False))
  for ix in range(nx):
    if is_good[iy,ix]>.5:
      # pull out basic parameters
      basicinfo = full_info[iy,ix,:pyirc.swi.Nb].tolist()
      #print('old current ->', basicinfo[pyirc.swi.I])
      basicinfo[pyirc.swi.I] = Ie[iy,ix]
      basicinfo[pyirc.swi.beta] = full_info[iy,ix,pyirc.swi.Nbb+1:pyirc.swi.Nbb+pyirc.swi.p] # in DN, +
      beta_in_e = -basicinfo[pyirc.swi.beta]/basicinfo[pyirc.swi.g]**numpy.linspace(1,pyirc.swi.p-1,num=pyirc.swi.p-1) # in e , -

      tslices0 = numpy.asarray([ts_vis, ts_vis+tchar1_vis, ts_vis+tchar2_vis])
      # initialize vector to stack correlation matrices:
      corr_stack = []
      for k in range(nvis):
        tslicesk = (tslices0+k).tolist()
        region_cube = pyirc.pixel_data(vislightfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], tslicesk,
                      [sensitivity_spread_cut, True], False)
        dark_cube = pyirc.pixel_data(visdarkfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], tslicesk,
                      [sensitivity_spread_cut, False], False)
        if fullref:
          lightref = lightref_array[k]
          darkref = darkref_array[k]
          #lightref = pyirc.ref_array(vislightfiles, formatpars, ny, tslicesk, False)
          #darkref = pyirc.ref_array(vislightfiles, formatpars, ny, tslicesk, False)
        else:
          lightref = numpy.zeros((len(vislightfiles), ny, 2*len(tslicesk)+1))
          darkref = numpy.zeros((len(visdarkfiles), ny, 2*len(tslicesk)+1))
        info = pyirc.corr_5x5(region_cube, dark_cube, tslicesk, lightref[:,iy,:], darkref[:,iy,:], basicpar, False)

        corr_matrix = info[4]
        var1 = info[2]
        var2 = info[3]
        # center of corr_matrix is element (2, 2) of the numpy array
        corr_matrix[2][2] = var2 - var1

        # median corrections to the central array of the auto-correlation matrix
        # (so we multiply the measured variance by the measured/predicted median,
        # this would perfectly correct for errors in Ie if the detector were exactly linear)
        med21 = info[1]
        predictmed = (tslicesk[2]*Ie[iy,ix]*(1. - numpy.sum(beta_in_e * (tslicesk[2]*Ie[iy,ix])**numpy.linspace(1,pyirc.swi.p-1,num=pyirc.swi.p-1)) )\
                     - tslicesk[1]*Ie[iy,ix]*(1. - numpy.sum(beta_in_e * (tslicesk[1]*Ie[iy,ix])**numpy.linspace(1,pyirc.swi.p-1,num=pyirc.swi.p-1)) ))\
                     / basicinfo[pyirc.swi.g]
        if basicpar.vis_med_correct: corr_matrix[2][2] /= med21/predictmed

        corr_stack.append(corr_matrix)
        # end loop over k

      corr_mean = numpy.mean(corr_stack, axis=0)
      # corr_mean is the v vector of eq. 34

      # now get the cube of data for BFE
      region_cube = pyirc.pixel_data(vislightfiles, formatpars, [dx*ix, dx*(ix+1), dy*iy, dy*(iy+1)], tslices,
                    [sensitivity_spread_cut, True], False)
    
      # iterate to solve BFE, Phi
    
      np2 = 2
      bfepar.Phi = numpy.zeros((2*np2+1,2*np2+1)); bfepar.Phi[np2,np2] = 1.e-12 # initialize to essentially zero
      if copy_ir_bfe:
        bfek_ir = full_info[iy,ix,pyirc.swi.Nb:pyirc.swi.Nbb].reshape((2*np2+1,2*np2+1))
        bfek = numpy.copy(bfek_ir)
      else:
        bfek  = pyirc.bfe(region_cube, tslices, basicinfo, bfepar, False) 
      tol = 1e-9
      diff = 1
      count = 0
      NN = 21
    
      while numpy.max(numpy.abs(diff)) > tol:   

        ts_vis_ref = ts_vis - basicpar.reset_frame
        tslices_vis = [ts_vis_ref,ts_vis_ref+tchar2_vis,ts_vis_ref,ts_vis_ref+tchar2_vis,nvis]
        tslices_vis1 = [ts_vis_ref,ts_vis_ref+tchar1_vis,ts_vis_ref,ts_vis_ref+tchar1_vis,nvis]
        normPhi = numpy.sum(bfepar.Phi) # this is omega/(1+omega)
        omega = normPhi / (1-normPhi)
        p2 = bfepar.Phi/normPhi
        sigma_a = 0.
        avals = [basicinfo[pyirc.swi.alphaV], basicinfo[pyirc.swi.alphaH], basicinfo[pyirc.swi.alphaD]] # (aV, aH, aD)
        truecorr = ftsolve.solve_corr_vis_many(bfek,NN,basicinfo[pyirc.swi.I],basicinfo[pyirc.swi.g],
                                       beta_in_e,sigma_a,tslices_vis,avals,omega=omega,p2=p2)
        #if count==0:
        #  print(tslices_vis, p2, truecorr)
        truecorr[2][2] = (truecorr-ftsolve.solve_corr_vis_many(bfek,NN,basicinfo[pyirc.swi.I],basicinfo[pyirc.swi.g],
                                       beta_in_e,sigma_a,tslices_vis1,avals,omega=omega,p2=p2))[2][2]
        diff = basicinfo[pyirc.swi.g]**2/(2*basicinfo[pyirc.swi.I]*tchar2_vis) * (corr_mean - truecorr)
        diff[2][2] = basicinfo[pyirc.swi.g]**2/(2*basicinfo[pyirc.swi.I]*(tchar2_vis-tchar1_vis)) * (corr_mean[2][2] - truecorr[2][2])
        bfepar.Phi += .5*(diff + numpy.flip(diff)) # force symmetrization here to avoid instability
    
        # update BFE
        if copy_ir_bfe:
          bfek = numpy.copy(bfek_ir)
        else:
          bfek  = pyirc.bfe(region_cube, tslices, basicinfo, bfepar, False) 
        count += 1
        
        if count>100:
            print('100 iterations of BFE/Phi solver reached, diff={:0.6f}'.format(numpy.max(numpy.abs(diff))))
            break

      #print('iter', count, 'omega = ',omega, 'max diff =', numpy.max(numpy.abs(diff)))

      # save information
      vis_bfek[iy,ix,:,:] = bfek
      vis_Phi[iy,ix,:,:] = bfepar.Phi
      op2 = ftsolve.op2_to_pars(bfepar.Phi)
      QYomega[iy,ix] = op2[0]
      cdCov[iy,ix,0] = op2[1]
      cdCov[iy,ix,1] = op2[2]
      cdCov[iy,ix,2] = op2[3]
      cdNiter[iy,ix] = op2[-1]

      # end loop over super-pixels
print('|')
print('')

# Now get ready to write information
print('Mean BFE kernel:')
print(numpy.mean(vis_bfek,axis=(0,1)))
print('Mean Phi kernel:')
print(numpy.mean(vis_Phi,axis=(0,1)))
print('sigma Phi kernel:')
print(numpy.std(vis_Phi,axis=(0,1)))
print('Charge diffusion parameters:')
print(ftsolve.op2_to_pars(numpy.mean(vis_Phi,axis=(0,1))))

# put all information into a gigantic array
vis_out_data = numpy.zeros((ny,nx,56))
vis_out_data[:,:,:25] = vis_bfek.reshape(ny,nx,25)
vis_out_data[:,:,25:50] = vis_Phi.reshape(ny,nx,25)
vis_out_data[:,:,50] = QYomega
vis_out_data[:,:,51:54] = cdCov
vis_out_data[:,:,54] = Ie
vis_out_data[:,:,55] = cdNiter
ncol = 56
#
# now we have in each super-pixel, 55 "columns" of data
# columns  0 .. 24 are the visible BFE kernel in e^-1 (order: dy=-2 dx=-2; dy=-2 dx=-1; dy=-2 dx=0; ...)
# columns 25 .. 49 are the visible Phi kernel (order: dy=-2 dx=-2; dy=-2 dx=-1; dy=-2 dx=0; ...)
# column 50 is the quantum yield omega parameter
# column 51 is Cxx charge diffusion in pixels^2
# column 52 is Cxy charge diffusion in pixels^2
# column 53 is Cyy charge diffusion in pixels^2
# column 54 is visible current Ie (e per frame)
# column 55 is number of iterations in p2 kernel

print ('')
print (vis_out_data.shape)
print ('Number of good regions =', numpy.sum(is_good))
mean_vis_out_data = numpy.mean(numpy.mean(vis_out_data, axis=0), axis=0)/numpy.mean(is_good)
std_vis_out_data = numpy.sqrt(numpy.mean(numpy.mean(vis_out_data**2, axis=0), axis=0)/numpy.mean(is_good) - mean_vis_out_data**2)
print('column, mean, stdev, stdev on the mean:')
for k in range(ncol):
  print('{:2d} {:12.5E} {:12.5E} {:12.5E}'.format(k, mean_vis_out_data[k], std_vis_out_data[k], std_vis_out_data[k]/numpy.sqrt(numpy.sum(is_good)-1)))
print ('')
#
# save to file
numpy.savetxt(outstem+'_visinfo.txt', vis_out_data.reshape(ny*nx, ncol))
numpy.save(outstem+'_visinfo.npy', vis_out_data)

# Saving some figures of these quantities:
matplotlib.rcParams.update({'font.size': 12})
num_bins = 30
F = plt.figure(figsize=(8,6))
S = F.add_subplot(2,2,1)
S.hist(QYomega.ravel(),bins=numpy.linspace(0, 0.1, num=num_bins))
S.set_xlabel(r'$\omega$')

S = F.add_subplot(2,2,2)
S.hist(Ie.ravel(),bins=num_bins)
S.set_xlabel(r'$I_e$')

S = F.add_subplot(2,2,3)
S.hist(cdNiter.ravel(),bins=numpy.linspace(0, 100, num=num_bins))
S.set_xlabel(r'Number of iterations')

S = F.add_subplot(2,2,4)
S.hist(vis_out_data[:,:,51].ravel(), num_bins, histtype='step', label=r'$C_{xx}$', linewidth=1.5, linestyle='-')
S.hist(vis_out_data[:,:,52].ravel(), num_bins, histtype='step', label=r'$C_{xy}$', linewidth=1.5, linestyle='--')
S.hist(vis_out_data[:,:,53].ravel(), num_bins, histtype='step', label=r'$C_{yy}$', linewidth=1.5, linestyle='-.')
S.set_xlabel(r'Charge diffusion component in pixels$^2$')
S.legend(loc='upper right', fontsize=12,frameon=False)
F.set_tight_layout(True)
F.savefig(outstem+'_vis_hist.pdf', bbox_inches='tight')
plt.close(F)

F = plt.figure(figsize=(15,8))
S = F.add_subplot(2,3,1)
S.set_title(r'$\omega$')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(QYomega, cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(2,3,2)
S.set_title(r'$I_e$')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
#S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(Ie, cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(2,3,3)
S.set_title(r'Number of iterations')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
#S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(cdNiter, cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(2,3,4)
S.set_title(r'$C_{xx}$')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,51], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(2,3,5)
S.set_title(r'$C_{xy}$')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
#S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,52], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(2,3,6)
S.set_title(r'$C_{yy}$')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
#S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,53], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

# F.set_tight_layout(True)
F.savefig(outstem+'_vis_matrices.pdf', bbox_inches='tight')
plt.close(F)

print('END')
