import numpy as np
import sys
from scipy.integrate import quad
import matplotlib.pyplot as plt

# get data
line = []
file1 = open(sys.argv[1], 'r') # minimizing_2.py qdat_F2_wavelengths.txt
count = 0
while True:
  count += 1
  l = file1.readline()
  if not l: break
  line += [l]
file1.close()

# define some constants
T = 95 #K
x = 0.445
Eg = -0.295 + 1.87*x - 0.28*x**2 + (6 - 14*x + 3*x**2)*(1E-4)*T + 0.35*x**4 #eV
alphag = (-65 + 1.88*T + (8694 - 10.31*T)*x)*100 #cm^-1 * 100 = m^-1
beta = -1 + 0.083*T + (21 - 0.13*T)*x 
hc = 1.23984193e-6 #eVm
alpha1 = .0167

# defined parameters that were found through minimization using 
h = 9.95877489e-07 #m 
xi = 1.00000000e+00 #unitless
params = [h, xi]

# create empty arrays & dictionary
wavelengths = np.array([])
my_dict = {}

# get experimental values
for l in line:
  arr = l.split()
  w = arr[0] # nm
  wavelengths = np.append(wavelengths, w) 
  for w in wavelengths:
      if w not in my_dict:
          my_dict[w] = [[],[],[]] # three lists: u, mtf2, mtf2err
  if len(arr)>4:
      if arr[0] in my_dict:
          u = float(arr[1]) # cyc/pix
          mtf2 = float(arr[-2])
          mtf2err = float(arr[-1])
          my_dict[arr[0]][0].append(u)
          my_dict[arr[0]][1].append(mtf2) 
          my_dict[arr[0]][2].append(mtf2err) 
  else:
      if arr[0] in my_dict:
          u = float(arr[1]) #cyc/pix
          mtf = float(arr[-2])
          mtferr = float(arr[-1])
          mtf2 = mtf**2
          mtf2err = mtf*2*mtferr
          my_dict[arr[0]][0].append(u)
          my_dict[arr[0]][1].append(mtf2)
          my_dict[arr[0]][2].append(mtf2err)
#print(my_dict)

def theory_mtf(z, params, u, w):
    # define the parameters
    h, xi = params
    # define some equations
    E = hc/w #eV
    alpha = (alphag * np.exp(np.sqrt(beta * (E - Eg)))) # m^-1
    # setting v = 0 
    r_p = 1/2 * (xi/h + np.sqrt((xi/h)**2 + 16*(np.pi)**2*u**2)) # m^-1
    r_m = 1/2 * (xi/h - np.sqrt((xi/h)**2 + 16*(np.pi)**2*u**2)) # m^-1
    P = (alpha * np.exp(-alpha*(h - z)))/(1 - np.exp(-1*alpha * h)) # m^-1
    G = ((((np.exp(r_p * z) - np.exp(r_m * z))*((r_m**2 * np.exp(r_m * z) * np.exp(h * (r_p - r_m))) - (r_p**2 * np.exp(r_p * z))))
          /(r_p * np.exp(r_p * z) - r_m * np.exp(r_m * z) * np.exp(h * (r_p - r_m)))) + (r_p * np.exp(r_p * z) - r_m * np.exp(r_m * z)))**(-1)*(r_p - r_m) 
    integrand = P*G
    return integrand

u_values = np.linspace(0, 0.98, 100) #cyc/pix
w = int(wavelengths[0])*1e-9 # 850 nm to m 
mtf2_values = []
#u_values = [u*1e5 for u in u_values] # converting u from cyc/pix to cyc/m
#print(u)
for u in u_values:
    mtf_theory, error = quad(theory_mtf, 0, h, args=(params, u*1e5, w), limit=5000)
    mtf2 = ((1-2*alpha1*(1-np.cos(2*np.pi*u)))*np.sinc(u)*mtf_theory)**2 # double scalars issue is coming from u=0, not a real concern
    mtf2_values = np.append(mtf2_values, mtf2)
    print(u, mtf_theory, mtf2)
    #print(mtf2)
#print(my_dict)     
# assigning the lists to plot
for key, value in my_dict.items():
    if key == '850':
        u_850 = value[0]
        mtf2_850 = value[1]
    if key == '980':
        u_980 = value[0]
        mtf2_980 = value[1]
    if key == '1310':
        u_1310 = value[0]
        mtf2_1310 = value[1]
    if key == '1550':
        u_1550 = value[0]
        mtf2_1550 = value[1]
    if key == '2000':
        u_2000 = value[0]
        mtf2_2000 = value[1]

plt.title("MTF determination", fontsize = 18) #y vs x
plt.xlabel('spatial frequency, u (cyc/pix)', fontsize = 12)
plt.ylabel('modulation transfer function, MTF$^2$', fontsize = 12)
plt.scatter(u_850, mtf2_850, c = 'gold')
plt.scatter(u_980, mtf2_980, c ='orange')
plt.scatter(u_1310, mtf2_1310, c='coral')
plt.scatter(u_1550, mtf2_1550, c='mediumvioletred')
plt.scatter(u_2000, mtf2_2000, c='darkviolet')
plt.plot(u_values, mtf2_values, c='darkblue')
plt.legend(['850nm', '980nm', '1310nm', '1550nm', '2000nm', 'robust model'], loc = 'upper right', fontsize = 12)
plt.savefig('MTF_theory')
plt.show()
    
