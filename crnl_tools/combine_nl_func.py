import numpy
import sys
import re
import matplotlib
import matplotlib.pyplot as plt

# maximum signal to use (in DN)
Spmax = 65536.

# main combination function
# takes in Smin[NF], Smax[NF], coefsIn[NF,ND+1]
# returns coefsOut[ND+1], a[NF], c[NF]
def combine_nl(Smin, Smax, coefsIn):
  NF = numpy.shape(coefsIn)[0]
  ND = numpy.shape(coefsIn)[1]-1

  Nv = NF*2+ND-1
  b = numpy.zeros((Nv))
  A = numpy.zeros((Nv,Nv))

  # useful later on
  x0 = numpy.poly1d([1])
  x1 = numpy.poly1d([1,0])

  for alpha in range(NF):
    f = numpy.poly1d(coefsIn[alpha,::-1]) #; print f
    dr = Smax[alpha] - Smin[alpha]

    # b coefficients
    #
    for j in range(2,ND+1):
      I = numpy.polyint(x1**(j+1))
      b[j-2] -= (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr
    #
    I = numpy.polyint(x1*f)
    b[ND-1 + alpha] += (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr
    #
    I = numpy.polyint(x1)
    b[ND-1 + NF + alpha] += (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr

    # A coefficients
    #
    for j in range(2,ND+1):
      for k in range(2,ND+1):
        I = numpy.polyint(x1**(j+k))
        A[j-2,k-2] += (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr
    #
    for j in range(2,ND+1):
      I = numpy.polyint(x1**j * f)
      A[j-2, ND-1 + alpha] -= (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr
      A[ND-1 + alpha, j-2] -= (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr
    #
    for j in range(2,ND+1):
      I = numpy.polyint(x1**j)
      A[j-2, ND-1 + NF + alpha] -= (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr
      A[ND-1 + NF + alpha, j-2] -= (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr
    #
    I = numpy.polyint(f**2)
    A[ND-1 + alpha, ND-1 + alpha] += (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr
    #
    I = numpy.polyint(f)
    A[ND-1 + alpha, ND-1 + NF + alpha] += (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr
    A[ND-1 + NF + alpha, ND-1 + alpha] += (numpy.polyval(I, Smax[alpha]) - numpy.polyval(I, Smin[alpha]))/dr
    #
    A[ND-1 + NF + alpha, ND-1 + NF + alpha] += 1.

  v = numpy.linalg.solve(A,b)
  outCoefs = numpy.zeros((ND+1))
  outCoefs[2:] = v[:ND-1]
  outCoefs[1] = 1.

  return outCoefs, v[ND-1:ND-1+NF], v[ND-1+NF:ND-1+2*NF]

# Now loop over channels

thisOut = open('nl.txt', 'w')
for channel in range(1,33):

  # Now read in the data
  files = [line.split(' ')[0]+'_offsets.txt' for line in open('currents.txt')]
  NumFiles = len(files)
  deg = 0
  for a in range(NumFiles):
    data = [line.rstrip('\n').lstrip(' ') for line in open(files[a])]
    usedata = re.split('\s+', data[channel-1])
    if a==0:
      deg = len(usedata)-2
      Sm = numpy.zeros((NumFiles))
      Sp = numpy.zeros((NumFiles))
      kja = numpy.zeros((NumFiles, deg+1))

    # fill in data
    Sm[a] = float(usedata[0])
    Sp[a] = float(usedata[1])
    if Sp[a]>Spmax: Sp[a]=Spmax
    kja[a,0] = 0.
    kja[a,1] = 1.
    for j in range(2,deg+1): kja[a,j] = float(usedata[1+j])

  print Sm; print ''; print Sp; print ''; print kja

  combined_poly, a_, c_ = combine_nl(Sm,Sp,kja)

  print 'Answer:', combined_poly
  print 'Offsets:', a_, c_

  for j in range(2,deg+1):
    thisOut.write(' {:19.12E}'.format(combined_poly[j]))
  thisOut.write('\n')

  # make plot
  e = .05
  matplotlib.rcParams.update({'font.size': 8})
  F = plt.figure(figsize=(6,8))
  S = F.add_subplot(2,1,1)
  S.set_title('NL curve: channel {:d}'.format(channel))
  S.set_xlabel('Signal [DN]')
  S.set_ylabel('Corrected - raw signal [DN]')
  S.set_xlim(-numpy.max(Sp)*e,numpy.max(Sp)*(1.+e))
  S.grid(True, color='g', linestyle='--')
  polya = numpy.poly1d(combined_poly[::-1])
  xa = numpy.linspace(numpy.min(Sm), numpy.max(Sp), 512)
  S.plot(xa,numpy.polyval(polya,xa)-xa,'r-')
  #
  S = F.add_subplot(2,1,2)
  S.set_title('NL curve residuals: channel {:d}'.format(channel))
  S.set_xlabel('Signal [DN]')
  S.set_ylabel('Stiching residual [linearized DN]')
  S.set_xlim(-numpy.max(Sp)*e,numpy.max(Sp)*(1.+e))
  S.grid(True, color='g', linestyle='--')
  #
  for a in range(NumFiles):
    xa = numpy.linspace(Sm[a],Sp[a],512)
    polya = numpy.poly1d(kja[a,::-1])
    polya *= a_[a]
    polya += c_[a]
    polya -= numpy.poly1d(combined_poly[::-1])
    S.plot(xa,numpy.polyval(polya,xa),'b-')
  #F.set_tight_layout(True)
  F.savefig(sys.argv[1]+'/nl-ch{:d}.pdf'.format(channel))
  plt.close(F)

  # end channel loop

thisOut.close

