import numpy as np
from ftsolve import solve_corr

def perturb_test(BFEK_model,diff,N,I,gain,beta,sigma_a,tslices,avals,avals_nl):
    dydx_max = np.zeros(BFEK_model.shape)
    for i in xrange(len(BFEK_model)):
	for j in xrange(len(BFEK_model[0])):
	    perturbed = np.copy(BFEK_model)
	    perturbed[i][j] += diff

	    solve_original = solve_corr(BFEK_model,N,I,gain,beta,sigma_a,tslices,avals,
		avals_nl)*((gain**2)/(I**2*(tslices[1]-tslices[0])
		*(tslices[-1]-tslices[-2])))

	    solve_perturbed = solve_corr(perturbed,N,I,gain,beta,sigma_a,tslices,avals,
		avals_nl)*((gain**2)/(I**2*(tslices[1]-tslices[0])
		*(tslices[-1]-tslices[-2])))

	    vals = ((solve_perturbed-solve_original)/diff)
  	    maxval = vals.flat[np.abs(vals).argmax()]
	    dydx_max[i][j] = maxval

    return dydx_max

