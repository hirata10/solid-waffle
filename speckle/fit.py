import numpy
import sys

# get data
line = []
file1 = open(sys.argv[1], 'r')
count = 0
while True:
  count += 1
  l = file1.readline()
  if not l: break
  line += [l]
file1.close()

alpha = .0167

sigma = numpy.linspace(0,.4,4001)
chi2 = numpy.zeros_like(sigma)
for l in line:
  arr = l.split()
  u = float(arr[0])
  if len(arr)>3:
    mtf2 = float(arr[-2])
    mtf2err = float(arr[-1])
  else:
    mtf = float(arr[-2])
    mtferr = float(arr[-1])
    mtf2 = mtf**2
    mtf2err = mtf*2*mtferr
  #print(u, mtf2, mtf2err)

  #theory_mtf = numpy.sinc(u) * (1-2*alpha*(1-numpy.cos(numpy.pi*u))) * numpy.exp(-2*numpy.pi**2*u**2*sigma**2)
  theory_mtf = numpy.sinc(u) * (1-2*alpha*(1-numpy.cos(numpy.pi*u))) / numpy.cosh(2*numpy.pi*u*sigma)
  chi2 += (mtf2-theory_mtf**2)**2/mtf2err**2

c2 = numpy.asarray(chi2)

for i in range(len(c2)): print(sigma[i], c2[i])

index = numpy.argmin(c2)
cdof = c2[index]/(len(line)-1)
c2 /= cdof
c2min = c2[index]
cLow = numpy.amin(numpy.where(c2<c2min+4)[0])
cHigh = numpy.amax(numpy.where(c2<c2min+4)[0])
print('{:6.4f} +{:6.4f}-{:6.4f}   {:6.3f}'.format(sigma[index], sigma[cHigh]-sigma[index], sigma[index]-sigma[cLow], cdof))

#print(chi2)
