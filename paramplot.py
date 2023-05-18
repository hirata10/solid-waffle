import numpy as np
import matplotlib
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
 
sim_gain = 2.06
sim_alpha = 1.69
sim_alphaH = 1.69
sim_alphaV = 1.69
sim_beta2 = -1.5725*1e-6
sim_beta3 = 1.9307e-5*1e-6
sim_beta4 = -1.4099e-10*1e-6
sim_ipnl = -1.1590
sim_ipnlNN = 0.2034

sim_truevals = [sim_alphaV,sim_alphaH,sim_beta2*sim_gain*-1e6,
		sim_beta3*sim_gain**2*-1e10,sim_beta4*sim_gain**3*-1e15,sim_gain,sim_ipnl,sim_ipnlNN]
 
 
def read_summary(fname):
    if fname is None:
        return None
        
    colnames = ['superX','superY','goodpix','raw_gain','gain_alpha','gain_alphabeta',
            'alphaH','alphaV','betaNL','q_per_t','alphaD','cH','cV','ipnl(-2,-2)',
            'ipnl(-1,-2)','ipnl(0,-2)','ipnl(1,-2)','ipnl(2,-2)','ipnl(-2,-1)',
            'ipnl(-1,-1)','ipnl(0,-1)','ipnl(1,-1)','ipnl(2,-1)','ipnl(-2,0)','ipnl(-1,0)',
            'ipnl(0,0)','ipnl(1,0)','ipnl(2,0)','ipnl(-2,1)','ipnl(-1,1)','ipnl(0,1)',
            'ipnl(1,1)','ipnl(2,1)','ipnl(-2,2)','ipnl(-1,2)','ipnl(0,2)','ipnl(1,2)',
            'ipnl(2,2)','t_intercept','beta2','beta3','beta4']

    table = pd.read_table(fname,delim_whitespace=True,index_col=(0,1),names=colnames,comment='#')

    # Adjust units as appropriate to this plot
    table['ipnlNN'] = 1e6*(table['ipnl(0,1)']+table['ipnl(0,-1)']
                    + table['ipnl(1,0)']+table['ipnl(-1,0)'])/4.0
    table['ipnl(0,0)'] *= 1e6
    table['alphaH'] *= 100
    table['alphaV'] *= 100
    table['alphaD'] *= 100

    table['beta2'] *= -1e6
    table['beta3'] *= -1e10
    table['beta4'] *= -1e15
    return table

font = {'family' : 'normal',
        'weight' : 'regular',
        'size'   : 13}

matplotlib.rc('font', **font)

prefix = '/users/PCON0003/cond0088/Projects/detectors/sw_outputs/PaperIII/'
files = [prefix+'chris_20663st_summary.txt',
	 prefix+'chris_20663st_128x16_summary.txt',
	 prefix+'chris_20663st-cub_summary.txt',
	 prefix+'chris_20663st-lo_summary.txt',
	 prefix+'chris_20663st-short_summary.txt',
         prefix+'chris_20663st-med_summary.txt',        
	 prefix+'chris_20828st_summary.txt',
	 prefix+'chris_20828st_128x16_summary.txt',
	 prefix+'chris_20828st-cub_summary.txt',
 	 prefix+'chris_20828st-lo_summary.txt',
	 prefix+'chris_20828st-short_summary.txt',
         prefix+'chris_20828st-med_summary.txt',
         prefix+'chris_20829st_summary.txt',
	 prefix+'chris_20829st_128x16_summary.txt',
	 prefix+'chris_20829st-cub_summary.txt',
	 prefix+'chris_20829st-lo_summary.txt',
	 prefix+'chris_20829st-short_summary.txt',
	 prefix+'chris_20829st-med_summary.txt',
	 prefix+'full_quart_nl_paperI_JG_pyircv25_16by16_summary.txt',
	 prefix+'full_quart_nl_paperI_JG_pyircv25_32by32_summary.txt']

# Set up figure
ylabels = ['SCA 20663, fiducial','SCA 20663, 128x16','SCA 20663, cubic CNL','SCA 20663, lo (1 3 4 6)','SCA 20663, short (5 7 8 10)','SCA 20663, med (3 6 7 10)',
	   'SCA 20828, fiducial','SCA 20828, 128x16','SCA 20828, cubic CNL','SCA 20828, lo (1 3 4 6)','SCA 20828, short (5 7 8 10)','SCA 20828, med (3 6 7 10)',
	   'SCA 20829, fiducial','SCA 20829, 128x16','SCA 20829, cubic CNL','SCA 20829, lo (1 3 4 6)','SCA 20829, short (5 7 8 10)','SCA 20829, med (3 6 7 10)']
sim_ylabels = ['simulations, 16x16','simulations, 32x32']
divisions = [5,11,17]
axis_names = [r'$\alpha_V$',r'$\alpha_H$',
              r'$\beta_2g$',r'$\beta_3g^2$',r'$\beta_4g^3$',
              r'$g$',r'$[K^2a+KK^I]_{0,0}$',r'$[K^2a+KK^I]_{<1,0>}$']
units = ['%','%',r'$10^6\times$DN$^{-1}$',r'$10^{10}\times$DN$^{-2}$',r'$10^{15}\times$DN$^{-3}$','e/DN',r'ppm/e',r'ppm/e']

xdata_labels = ['alphaV','alphaH','beta2','beta3','beta4',
                'gain_alphabeta','ipnl(0,0)','ipnlNN']

colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:(divisions[1]-divisions[0])]

fsize = (2*len(axis_names),15)

fig = plt.figure(figsize=fsize)
grid = plt.GridSpec(15,len(axis_names))
axes=[]
sim_axes = []
for i in range(len(axis_names)):
    ax = fig.add_subplot(grid[:-2,i])
    axes.append(ax)

    sax = fig.add_subplot(grid[-2:,i])
    sim_axes.append(sax)

tables = [read_summary(f) for f in files]

for i,ax in enumerate(axes):
    ax.set_title(axis_names[i])
    ax.set_ylim(len(ylabels)-0.5,-0.5)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(-3,3))
    ax.set_yticklabels([])
    
    if i==0:
        ax.set_yticks([n for n in range(len(ylabels))])
        ax.set_yticklabels(ylabels)
        ax.set_ylim(len(ylabels)-0.5,-0.5)

    ax.tick_params(axis='y',which='both',left=False)
    ax.set_xlabel(units[i])    
    ax.xaxis.labelpad=15

    for j,run in enumerate(ylabels):
        mean = np.mean(tables[j][xdata_labels[i]])
        stdev = np.std(tables[j][xdata_labels[i]])
	error = stdev#/len(tables[j]**0.5)
        ax.errorbar([mean],[j],xerr=[error],capsize=8.0,markersize=10,marker='o',
                    color=colors[j%(divisions[1]-divisions[0])])
        if j==0:
            for k, div in enumerate(divisions):
                ax.axhline(y=div+0.5,color='k',linewidth=0.75,linestyle='-')

        if j%(divisions[1]-divisions[0])==0:
            ax.fill_betweenx([j-0.5,j+divisions[1]-divisions[0]-0.5],[mean-error,mean-error],
                             [mean+error,mean+error], color = 'grey', alpha = 0.5)
            
sim_tables = tables[-1*len(sim_ylabels):]

for i,sax in enumerate(sim_axes):
    sax.set_ylim(len(sim_ylabels)-0.5,-0.5)
    sax.ticklabel_format(axis='x', style='sci', scilimits=(-3,3))
    sax.set_yticklabels([])
    
    if i==0:
        sax.set_yticks([n for n in range(len(sim_ylabels))])
        sax.set_yticklabels(sim_ylabels)
        sax.set_ylim(len(sim_ylabels)-0.5,-0.5)

    sax.tick_params(axis='y',which='both',left=False)
    
    for j,run in enumerate(sim_ylabels):
	mean = np.mean(sim_tables[j][xdata_labels[i]])
	stdev = np.std(sim_tables[j][xdata_labels[i]])
        error = stdev#/len(sim_tables[j]**0.5)
	sax.errorbar([mean],[j],xerr=[error],capsize=8.0,markersize=10,marker='o',
                    color=colors[j%(divisions[1]-divisions[0])])

        if j%(divisions[1]-divisions[0])==0:
	    sax.axvline(x=sim_truevals[i],linestyle='--',color='k')
            #sax.fill_betweenx([j-0.5,j+len(sim_ylabels)-0.5],[mean-stdev,mean-stdev],
            #                 [mean+stdev,mean+stdev], color = 'grey', alpha = 0.5)

plt.tight_layout()

plt.savefig('paramplot_test.pdf',format='pdf')
plt.show()

