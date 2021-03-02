import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import SymLogNorm
from SCA_data_table import read_summary

fname = '/users/PCON0003/cond0088/Projects/detectors/sw_outputs/chris_20829st_summary.txt'

table = read_summary(fname)

rows = ['q_per_t','gain_alphabeta','alpha','alphaH','alphaV','alphaD',
        'beta2g','beta3g2','beta4g3','ipnl(0,0)','ipnl<1,0>','ipnl<1,1>',
        'ipnl<2,0>','ipnl<2,1>','ipnl<2,2>','ipnlH','ipnlV']

labels = [r'$It_{n,n+1}$',r'$g$',r'$\alpha$',r'$\alpha_\mathrm{H}$',
          r'$\alpha_\mathrm{V}$',r'$\alpha_\mathrm{D}$',r'$\beta_2 g$',
	  r'$\beta_3 g^2$',r'$\beta_4 g^3$',r'$[K^2a+KK^I]_{0,0}$',
	  r'$[K^2a+KK^I]_{\langle 1,0 \rangle}$',r'$[K^2a+KK^I]_{\langle 1,1 \rangle}$',
	  r'$[K^2a+KK^I]_{\langle 2,0 \rangle}$',r'$[K^2a+KK^I]_{\langle 2,1 \rangle}$',
	  r'$[K^2a+KK^I]_{\langle 2,2 \rangle}$',r'$[K^2a+KK^I]_\mathrm{H}$',
	  r'$[K^2a+KK^I]_\mathrm{V}$']

rcParams['figure.figsize'] = 10,10

im = plt.matshow(table.loc[:,rows].corr(),cmap='RdBu')
plt.xticks(rotation=45)
plt.xticks(range(len(rows)),labels,fontsize=14,rotation=90)
plt.yticks(range(len(rows)),labels,fontsize=14)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=14)

plt.savefig('param_corr.pdf',format='pdf')
plt.show()
