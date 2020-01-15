import numpy as np
from numpy.fft import fft2,ifft2
import warnings
import pdb

def center(arr):
    # Transforms kernel so that it looks as expected to the eye:
    # Returns version of kernel with (0,0) in center, y=-1 down, x=-1 left, etc.
    # Centered arrays should *NOT* be used for calculations! For display only :)
    size = (len(arr)+1)//2
    return np.roll(np.roll(arr,-size,axis=0),-size,axis=1)

def decenter(arr):
    # Transforms kernel from human-readable to numpy-readable 
    # Only decentered arrays should be used for calculations
    size = (len(arr)+1)//2
    return np.roll(np.roll(arr,size,axis=0),size,axis=1)

def flip(arr):
    # transforms decentered a_i,j -> decentered a_-i,-j
    arrc = center(arr)
    arrc_flipped = np.flip(arrc.flatten(),axis=0).reshape(arrc.shape)
    return decenter(arrc_flipped)

def pad_to_N(arr,N):
    # pads array out to size NxN (if arr is smaller than this)
    # input assumed to be centered
    if not arr.shape[0]>N:
        pad_size = (N-arr.shape[0])//2
        return np.pad(arr,pad_size,mode='constant')
    else:
        return arr

def solve_corr(bfek,N,I,g,betas,sigma_a,tslices,avals,avals_nl=[0,0,0],outsize=2):
    # INPUT: 
    # bfek     <- compound kernel [K^2 a+KK*](assumed to be centered)
    # N        <- detector size (assumed odd for now)
    # I        <- current
    # g        <- gain (assuming no higher order fitting for now)
    # betas    <- array of classical non-linearity coefficients [beta_2...beta_n]
    # sigma_a  <- sum of the BFE kernel
    # tslices  <- list of time slices (ta, tb, tc, td)  
    # avals    <- list of alpha values for linear IPC kernel (aV, aH, aD)
    # avals_nl <- list of alpha values for NL-IPC kernel (aV_nl, aH_nl, aD_nl)
    # outsize  <- "radius" of output (equiv. to sBFE_out in main code)
    # OUTPUT: C_abcd
    
    ta, tb, tc, td = tslices
    aV, aH, aD = avals
    aV_nl, aH_nl, aD_nl = avals_nl
    
    # convert betas to an array if it isn't already
    if not isinstance(betas, np.ndarray): betas = np.array([betas])

    if not bfek.shape[1]==bfek.shape[0]:
        warnings.warn("WARNING: convolved BFE kernel (BFEK) not square.")

    assert(N==2*(N//2)+1)
    
    # Calculate K and K* from given alphas
    cent = slice(N//2-outsize,N//2+outsize+1)

    k = decenter(pad_to_N(np.array([[aD,aV,aD],
              [aH,1-4*aD-2*aV-2*aH,aH],
              [aD,aV,aD]]),N))

    knl = decenter(pad_to_N(np.array([[aD_nl,aV_nl,aD_nl],
                [aH_nl,-4*aD_nl-2*aV_nl-2*aH_nl,aH_nl],
                [aD_nl,aV_nl,aD_nl]]),N))

    # solve Fourier version for asq: F(BFEK) = Ksq^2*asq + Ksq*Knl_sq  
    bfek = decenter(pad_to_N(bfek,N))
    ksq = fft2(k)
    knl_sq = fft2(knl)
    asq = (fft2(bfek)- ksq*knl_sq)/ksq**2
    a = ifft2(asq)

    a_flipped = flip(a)
    afsq = fft2(a_flipped)
    afsq_p = flip(afsq)

    ksq_p = flip(ksq)
    knl_sq_p = flip(knl_sq)

   # Calculate Cov(qsq(t),qsq(t')) (see eqn 38)
    qqs = []

    for ts in [(ta,tc),(ta,td),(tb,tc),(tb,td)]:
        t1 = min(ts)
        t = max(ts)
    
        #qq = (1/(afsq+afsq_p-sigma_a) * np.exp(I*afsq*(t-t1)) *
         #   (np.exp(I*(afsq+afsq_p)*t1)-np.exp(I*sigma_a*t1)))
        
        X = I*t1*(afsq+afsq_p-sigma_a)
        qq = (np.where(np.abs(X)>1e-4, (np.exp(X)-1)/np.where(np.abs(X)>1e-5,X,X+1),
                          1+X/2.+X**2/6.+X**3/24.))*I*t1*np.exp(I*afsq*(t-t1))*np.exp(I*sigma_a*t1)                          

        qqs.append(qq)
    
# Plug into correlation function (see eqn 51)
    csq_abcd =(1/g**2
           *(eval_cnl(betas,I,ta)*eval_cnl(betas,I,tc)*(ksq+knl_sq*I*ta)*(ksq_p+knl_sq_p*I*tc)*qqs[0] 
           - eval_cnl(betas,I,ta)*eval_cnl(betas,I,td)*(ksq+knl_sq*I*ta)*(ksq_p+knl_sq_p*I*td)*qqs[1] 
           - eval_cnl(betas,I,tb)*eval_cnl(betas,I,tc)*(ksq+knl_sq*I*tb)*(ksq_p+knl_sq_p*I*tc)*qqs[2] 
           + eval_cnl(betas,I,tb)*eval_cnl(betas,I,td)*(ksq+knl_sq*I*tb)*(ksq_p+knl_sq_p*I*td)*qqs[3])
           )
    
    return center(np.real(ifft2(csq_abcd)))[cent][:,cent]

def eval_cnl(betas,I,t):
    nu = np.arange(2,len(betas)+2)
    return 1-np.sum(nu*betas*(I*t)**(nu-1))


# same as solve_corr *except* that we have tslice [ta,tb,tc,td,tn],
# where tn >= 1 is the number of similar steps to use -- i.e., we have
# (C_{ta,tb,tc,td} + C_{ta+1,tb+1,tc+1,td+1} + ... + C_{ta+tn-1,tb+tn-1,tc+tn-1,td+tn-1} )/tn
def solve_corr_many(bfek,N,I,g,betas,sigma_a,tslices,avals,avals_nl=[0,0,0],outsize=2):
   this_t = tslices[:-1]
   tn = tslices[-1]
   cf = solve_corr(bfek,N,I,g,betas,sigma_a,this_t,avals,avals_nl,outsize)
   for j in range(tn-1):
     for k in range(4): this_t[k] += 1
     cf += solve_corr(bfek,N,I,g,betas,sigma_a,this_t,avals,avals_nl,outsize)
   cf /= tn+0.0
   return(cf)

#if __name__=="__main__":
   
   # Test against configuration-space corrfn generated from known inputs/simulated flats
#   N = 21
#   I = 1487
#   g = 2.06
#   beta = 0
#   tslices = [3, 11, 13, 21]
#   avals = [0,0,0]
#   avals_nl = [0,0,0]   

#   test_bfek = 1.E-6*np.array(
#    	[[-0.01, 0.0020, -0.0210, -0.019, 0.028],
#     	[0.0040, 0.0490, 0.2480, 0.01, -0.0240],
#     	[-0.0170, 0.2990, -1.372, 0.2840, 0.0150],
#     	[0.0130, 0.0560, 0.2890, 0.0390, 0.02],
#     	[0.035, 0.0070, 0.0380, 0.0010, 0.026]])


   #test_bfek = np.load('/users/PCON0003/cond0088/Projects/detectors/solid-waffle/testBFEK_flatsim_matcheddark_bfeonly18237sim_10files_sub20.npy')
#   sigma_a = np.sum(test_bfek)

   # Test against BFEK values in run of test_run.py with input config.18237.sample1
   #N = 21
   #I = 1378
   #g = 2.26
   #beta = 5.98e-7
   #sigma_a = 0.0
   #tslices = [3, 11, 13, 21]
   #avals = [0.014,0.023,0]
   #avals_nl = [0,0,0]
   #test_bfek = np.load('test_bfek.npy')

#   c_abcd = solve_corr(test_bfek,N,I,g,beta,sigma_a,tslices,avals,avals_nl)
#   print (c_abcd)

