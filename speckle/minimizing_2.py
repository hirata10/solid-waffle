import numpy as np
import sys
from scipy.integrate import quad
from scipy.optimize import minimize

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
alphag = (-65 + 1.88*T + (8694 - 10.31*T)*x)*100 #m^-1
beta = -1 + 0.083*T + (21 - 0.13*T)*x
hc = 1.23984193e-6 #eVm
alpha1 = .0167

# parameters to be minimized over - best guesses
h = 1e-06 #m
xi = 1 #unitless
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
          u = float(arr[1]) #cyc/pix
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
    
# define the theory function for the mtf
def theory_mtf(z, params, u, w):
    # define the parameters
    h, xi = params
    # define some equations
    E = hc/w #eV
    alpha = alphag * np.exp(np.sqrt(beta * (E - Eg))) #m^-1
    # setting v = 0 
    r_p = 1/2 * (xi/h + np.sqrt((xi/h)**2 + 16*(np.pi)**2*u**2)) #m^-1
    r_m = 1/2 * (xi/h - np.sqrt((xi/h)**2 + 16*(np.pi)**2*u**2)) #m^-1
    P = (alpha * np.exp(-alpha*(h - z)))/(1 - np.exp(-1*alpha * h)) #m^-1
    G = ((((np.exp(r_p * z) - np.exp(r_m * z))*((r_m**2 * np.exp(r_m * z) * np.exp(h * (r_p - r_m))) - (r_p**2 * np.exp(r_p * z))))
          /(r_p * np.exp(r_p * z) - r_m * np.exp(r_m * z) * np.exp(h * (r_p - r_m)))) + (r_p * np.exp(r_p * z) - r_m * np.exp(r_m * z)))**(-1)*(r_p - r_m)   #unitless
    integrand = P*G #m^-1
    return integrand

def chi_squared(params, my_dict, theory_mtf):
    chi2 = np.array([])
    params_list = params.tolist()
    for key, value in my_dict.items():
        w = int(key)*1e-9 #m
        u_values = value[0] #cyc/pix
        mtf2_values = value[1]
        mtf2err_values = value[2]
        for u, mtf2, mtf2err in zip(u_values, mtf2_values, mtf2err_values):
            mtf_theory, error = quad(theory_mtf, 0, h, args=(params_list, u*1e5, w), limit=5000) # u is now cyc/m
            mtf2_theory = (np.sinc(u)*(1-2*alpha1*(1-np.cos(2*np.pi*u)))*mtf_theory)**2
            chi2 = np.append(chi2, (mtf2-mtf2_theory**2)**2/mtf2err**2)
    chi_squared_value = np.sum(chi2)
    return chi_squared_value

for tier in range(2):
    result = minimize(chi_squared, params, args=(my_dict, theory_mtf), options = {'disp':True, 'return_all':True,'scale': np.array([1e-6,1]),'maxfun': 1000 ,'initial_simplex': np.array([[6e-6,-.75], [5e-6,-.75], [6e-6,-1.]])}, method = 'TNC') #Nelder-Mead, TNC
    optimized_params = result.x
    print(optimized_params)
    print(chi_squared(optimized_params, my_dict, theory_mtf))
    params = optimized_params
