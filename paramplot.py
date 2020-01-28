import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
         
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
    table['ipnlNN'] = (table['ipnl(0,1)']+table['ipnl(0,-1)']
                    + table['ipnl(1,0)']+table['ipnl(-1,0)'])/4.0
    return table

prefix = '/users/PCON0003/cond0088/Projects/detectors/sw_outputs/'
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
	 prefix+'chris_20828st_128x16_summary.txt',
	 prefix+'chris_20829st-cub_summary.txt',
	 prefix+'chris_20829st-lo_summary.txt',
	 prefix+'chris_20829st-short_summary.txt',
	 prefix+'chris_20829st-med_summary.txt',
	 prefix+'full_quart_nl_paperI_JG_pyircv25_16by16_summary.txt',
	 prefix+'full_quart_nl_paperI_JG_pyircv25_32by32_summary.txt']

# Set up figure
ylabels = ['SCA 20663, fiducial','SCA 20663, 128x16','SCA 20663, cubic CNL','SCA 20663, lo (1 3 4 6)','SCA 20663, short (5 7 8 10)','SCA 20663, med (3 6 7 10)',
	   'SCA 20828, fiducial','SCA 20828, 128x16','SCA 20828, cubic CNL','SCA 20828, lo (1 3 4 6)','SCA 20828, short (5 7 8 10)','SCA 20828, med (3 6 7 10)',
	   'SCA 20829, fiducial','SCA 20829, 128x16','SCA 20829, cubic CNL','SCA 20829, lo (1 3 4 6)','SCA 20829, short (5 7 8 10)','SCA 20829, med (3 6 7 10)',
	   'simulations, 16x16','simulations, 32x32']
divisions = [5,11,17,23]
axis_names = [r'$\alpha_V$',r'$\alpha_H$',
              r'$\beta_2$',r'$\beta_3$',r'$\beta_4$',
              r'$g$',r'$[K^2a+KK\']_{0,0}$',r'$[K^2a+KK\']_{NN}$']
units = ['%','%',r'DN$^{-1}$',r'DN$^{-2}$',r'DN$^{-3}$','e/DN',r'e$^{-1}$',r'e$^{-1}$']

xdata_labels = ['alphaV','alphaH','beta2','beta3','beta4',
                'gain_alphabeta','ipnl(0,0)','ipnlNN']

colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:(divisions[1]-divisions[0])]
xlims = []
fsize = (3*len(axis_names),10)

fig,axes = plt.subplots(1,len(axis_names),figsize=fsize,sharey=True)

#fig.subplots_adjust(wspace=0)

for i,ax in enumerate(axes):
    ax.set_title(axis_names[i])
    ax.set_ylim(len(ylabels)-0.5,-0.5)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(-3,3))
    
    if i==0:
        #ax.spines['right'].set_visible(False)
        ax.set_yticks([i for i in range(len(ylabels))])
        ax.set_yticklabels(ylabels)
   # elif i==len(axes)-1:
        #ax.spines['left'].set_visible(False)
    #else:
        #ax.spines['right'].set_visible(False)
        #ax.spines['left'].set_visible(False)

    ax.tick_params(axis='y',which='both',left=False)
    
# Read data
tables = [read_summary(f) for f in files]

# Plot data
for i,ax in enumerate(axes):
    for j,run in enumerate(ylabels):
        mean = np.mean(tables[j][xdata_labels[i]])
        stdev = np.std(tables[j][xdata_labels[i]])
        ax.errorbar([mean],[j],xerr=[stdev],capsize=8.0,markersize=10,marker='o',
                    color=colors[j%(divisions[1]-divisions[0])])
        ax.set_xlabel(units[i])
	ax.xaxis.labelpad = 15

	if j==0:
            for k, div in enumerate(divisions):
		ax.axhline(y=div+0.5,color='k',linewidth=0.75,linestyle='-')		
	
	if j%(divisions[1]-divisions[0])==0:
	    ax.fill_betweenx([j-0.5,j+divisions[1]-divisions[0]-0.5],[mean-stdev,mean-stdev],
			     [mean+stdev,mean+stdev], color = 'grey', alpha = 0.7)

plt.savefig('paramplot.pdf',format='pdf')       
plt.show()


