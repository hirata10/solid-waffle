""" Functions to create various detector effect.
The structure will change, but for now this will be the location for 
functions related to IPC, BFE, etc.
"""
import numpy as np
from numpy.random import randn,poisson
sys.path.insert(0, '../')
#sys.path.insert(0, '/users/PCON0003/cond0080/src/solid-waffle/')
from pyirc import *

def simple_ipc_kernel(alpha=0.01):
  """ Simple function to return a 3 x 3 kernel with an alpha where alpha
  is the kernel value for the 4 adjacent pixels, and the central value is
  1-4*alpha.  This is a symmetric kernel.
  alpha typically on a percent level.
  """
  kernel = np.zeros((3, 3))
  kernel[1,0] = kernel[0,1] = kernel[1,2] = kernel[2,1] = alpha
  kernel[1,1] = 1-4*alpha
  return kernel
