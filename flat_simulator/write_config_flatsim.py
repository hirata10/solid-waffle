"""
script to write a config file for multiple noise realizations on the
commandline or in a bash script

hard-coded items include script and output locations and settings for
all the other parameters besides whether you want a flat or a dark and
noise realization

the choices of params depend on the detector you want to simulate, so one
should match things like gain, quantum efficiency, etc. for that detector
"""
import os
import glob
import datetime
import numpy as np
import sys

# User inputs
if len(sys.argv) != 5:
  print("Didn't give arguments to specify flat, name, rng")
  print("First arg specify 'flat', or any other arg if you want a dark")
  print("Second arg is a string to use in the names of the configs and sims")
  print("Third arg is a random number for the noise realization")
  print("Fourth arg specify 'ir' or 'vis' --- latter enables quantum yield")
  print("Example: python write_config_flatsim.py x test 10 ir")
  exit()
else:
  flatflag = sys.argv[1]
  name = str(sys.argv[2])
  rng = sys.argv[3]
  waveflag = sys.argv[4]

# Inits hard-coded to directories at OSC
tdir = '/users/PCON0003/cond0080/src/solid-waffle/flat_simulator/'
outdir = '/fs/project/PCON0003/ami/simulated_detector/'
today = datetime.date.today()
# Write out the first config file
if flatflag=='flat':
  file = open('%s/%ssim_config_%s_rng%d' %(tdir,waveflag,name,int(rng)),'w')
else:
  file = open('%s/%ssim_config_%s_%s_rng%d' %(
      tdir,waveflag,name,flatflag,int(rng)),'w')
    
file.write('# Generated %s\n' %str(today))
file.write('# Format (1 = H4RG, WFIRST-like, 2 is an H2RG)\n') 
file.write('FORMAT: 1\n')
file.write('NREADS: 64\n')  # Earlier prototype SCAs had 66 reads
file.write('SUBSTEPS: 20\n')
file.write('DT: 2.75\n')
file.write('GAIN: 1.7285\n')  # Gain similar to SCA 20829
if waveflag=='ir':
  pass
elif waveflag=='vis':
  # Format is QY_omega sig and then ftsolve.p2kernel takes as first arg cov
  # list or np array, length 3: xx, xy, yy
  # numbers hard-coded to what was used for Givans+ charge diffusion paper
  file.write('QY: 0.08 0.04 0.02 0.04\n')
else:
  sys.exit("ir vs vis was not specified, dying")
if flatflag=='flat':
  file.write('ILLUMINATION: 1109.64\n')
else:
  file.write('ILLUMINATION: 0.08\n')
file.write('QE: 9.5e-1\n')
file.write('RNGSEED: %d\n'%int(rng))
file.write('RESET_E: 2.14e4\n')
file.write('NOISE: Gauss\n')
file.write('WAVEMODE: %s\n'%str(waveflag))
file.write('BFE: true\n')
# Uncomment below if linear IPC is being set
#file.write('L_IPC: true\n')
# Uncomment below if non-linearity is being set
#file.write('NL: quartic -1.5725 1.9307e-5 -1.4099e-10\n')
if flatflag=='flat':
  file.write('OUTPUT: %s/flat_%s_%s_rng%d.fits\n'%(
      outdir,waveflag,name,int(rng)))
else:
  file.write('OUTPUT: %s/%s_%s_%s_rng%d.fits\n'%(
      outdir,flatflag,waveflag,name,int(rng)))

file.close()
