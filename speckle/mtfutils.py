import numpy

# functions to fit two triangles to a power spectrum
# (array, 0 element corresponding to 0 wavenumber)
#
# as a function of offset (cpp0) and width (cppw)
#
# output = best amplitude1, best amplitude 2, SSR
def get_ps_err(ps_in, cpp0, cppw):
  cpp0 += 1.
  n = numpy.size(ps_in)
  ps_inS = ps_in - numpy.median(ps_in)
  template1 = numpy.zeros((n,))
  template2 = numpy.zeros((n,))
  dx = cppw*n
  for j in range(1+int(dx)): template1[j] = template1[-j] = 1.-j/dx
  xc = cpp0*n
  xmin = xc-dx
  xmax = xc+dx
  for j in range(int(xmin),int(xmax)+1):
    amp = .5*(1.-numpy.abs(j-xc)/dx)
    template2[j%n] += amp
    template2[-(j%n)] += amp
  F11 = numpy.sum(template1[1:]**2)
  F12 = numpy.sum(template1[1:]*template2[1:])
  F22 = numpy.sum(template2[1:]**2)
  S1 = numpy.sum(ps_inS[1:]*template1[1:])
  S2 = numpy.sum(ps_inS[1:]*template2[1:])
  bestamp1 = (F22*S1-F12*S2)/(F11*F22-F12**2)
  bestamp2 = (F11*S2-F12*S1)/(F11*F22-F12**2)
  ssr = numpy.sum((ps_inS[1:]-bestamp1*template1[1:]-bestamp2*template2[1:])**2)
  return bestamp1,bestamp2,ssr

# return "best" power spectrum from triangel fit
# and normalized residuals:
#
# cpp0 best, cppw best, amplitude 1, amplitude 2, SSR/SS
def get_triangle_from_ps(ps_in, cppguess):

  gr = 33
  sp = 8

  nn = numpy.size(ps_in)
  ssr_array = numpy.zeros((gr,gr))
  cpp0 = numpy.linspace(1.05/nn,.5,gr)
  if cppguess>.4: cpp0 = numpy.linspace(.2,.5,gr)
  if cppguess>.5: cpp0 = numpy.linspace(.5,.8,gr)
  cppw = numpy.linspace(1.05/nn,.255,gr)
  for iter in range(24):
    for i in range(gr):
      for j in range(gr):
        a1,a2,ssr = get_ps_err(ps_in, cpp0[i], cppw[j])
        ssr_array[i,j] = ssr
    ind = numpy.unravel_index(numpy.argmin(ssr_array), numpy.shape(ssr_array))
    ii = ind[0]; jj = ind[1]
    best_cpp0 = cpp0[ii]; best_cppw = cppw[jj]

    # make next grid
    if ii<sp: ii=sp
    if jj<sp: jj=sp
    if ii>=gr-sp: ii=gr-sp-1
    if jj>=gr-sp: jj=gr-sp-1
    cpp0 = numpy.linspace(cpp0[ii-sp],cpp0[ii+sp],gr)
    cppw = numpy.linspace(cppw[jj-sp],cppw[jj+sp],gr)

  # flip if guess on other side
  a1, a2, ssr = get_ps_err(ps_in, best_cpp0, best_cppw)

  return best_cpp0, best_cppw, a1, a2, ssr/numpy.sum(ps_in[1:]**2)

# smothed triangle function:
# length L (x=0 ... x=L-1)
# peak at x=x0
# half width dx0
# smooth with |sin(pi x')/sin(pi x'/N)|^2/N
#
# ovsamp = oversampling factor
def trismooth(N, x0, dx0, ovsamp=16):

  NN = N*ovsamp
  x_ = numpy.linspace(0, (1-1/NN)*N, NN)
  S1 = (numpy.sinc(x_)/numpy.sinc(x_/N))**2*N
  xpk = x0 - N*numpy.floor(x0/N)
  S2 = numpy.maximum(0, 1-numpy.abs((x_-xpk)/dx0)) + numpy.maximum(0, 1-numpy.abs((x_-xpk-N)/dx0))
  SF = numpy.real(numpy.fft.ifft(numpy.fft.fft(S1)*numpy.fft.fft(S2)))/NN
  return SF[::ovsamp]

# create the phi matrix using characteristics of the power spectrum
#
# subregion -> only fits the region around the expected peak
# smoothtri -> uses sinc^2-smoothed instead of sharp triangles
#
# return the phi matrix
def make_phi_matrix(PS,u,delta_u,subregion=False,smoothtri=True):

  # make the empty array, with the length of the power spectra
  phi = numpy.empty((0,len(PS)))
  i_arr = numpy.array(range(len(PS))).astype(numpy.float64)

  # convert u to array values
  i = u * len(PS)
  delta_i = delta_u * len(PS)
  if not smoothtri: delta_i = numpy.sqrt(delta_i**2 + .886**2) # RSS addition of sinc^2 to FWHM
  offset = int(int(i) - len(PS)/2)

  # create the first triangle
  A = 1 - numpy.abs(i_arr - (i - offset))/delta_i
  A = numpy.maximum(A, 0)
  A = numpy.roll(A, offset)
  one = numpy.roll(A, int(len(PS) - i))

  one = numpy.maximum(numpy.roll(1 - numpy.abs(i_arr - len(PS)/2)/delta_i, len(PS)//2), 0.)

  if smoothtri:
    one = trismooth(len(PS),0.,delta_i)

  # add to the matrix
  phi = numpy.vstack((phi, one))

  # create the second triangle
  B = 1 - numpy.abs(i_arr - (len(PS) - (i - offset)))/delta_i
  B = numpy.maximum(B, 0)
  B = numpy.roll(B, len(PS) - offset)
  two = A + B

  if smoothtri:
    two = trismooth(len(PS),i,delta_i) + trismooth(len(PS),len(PS)-i,delta_i)

  # add to the matrix
  phi = numpy.vstack((phi, two))

  # integral constraint
  C = numpy.zeros(len(PS))
  C[0] = 1

  # add to the matrix
  phi = numpy.vstack((phi, C))

  # poisson noise
  P = numpy.zeros(len(PS))
  P[:] = 1

  # add to the matrix
  phi = numpy.vstack((phi, P))

  if subregion:
    i_ = len(PS) * numpy.modf(i/len(PS)+4)[0]
    im = min(i_,len(PS)-i_)
    ip = max(i_,len(PS)-i_)
    imin = int(numpy.floor(im-3*(delta_i+1)))
    imax = int(numpy.ceil(im+3*(delta_i+1)))
    if imax>len(PS)//2: imax = int(numpy.ceil(ip+3*(delta_i+1)))
    imin = max(imin,0)
    imax = min(imax,len(PS)-1)
    phi = two
    for j in range(3):
      phi = numpy.vstack((phi, ((i_arr-i_)/delta_i)**j))
    if imin<delta_i or smoothtri: phi = numpy.vstack((phi, one))
    # now symmetrize the whole thing
    phi[:,:imin] = 0
    phi[:,imax+1:] = 0
    # eliminate zero
    phi[:,0] = 0.

  return phi
