import numpy as np
import matplotlib
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16

def read_summary(fname):
    if fname is None:
        return None
        
    colnames = ['superX','superY','goodpix','raw_gain','gain_alpha','gain_alphabeta',
                'alphaH','alphaV','betaNL','q_per_t','alphaD','cH','cV',
                'ipnl(-2,-2)','ipnl(-1,-2)','ipnl(0,-2)','ipnl(1,-2)','ipnl(2,-2)',
                'ipnl(-2,-1)', 'ipnl(-1,-1)','ipnl(0,-1)','ipnl(1,-1)','ipnl(2,-1)',
                'ipnl(-2,0)','ipnl(-1,0)','ipnl(0,0)','ipnl(1,0)','ipnl(2,0)',
                'ipnl(-2,1)','ipnl(-1,1)','ipnl(0,1)','ipnl(1,1)','ipnl(2,1)',
                'ipnl(-2,2)','ipnl(-1,2)','ipnl(0,2)','ipnl(1,2)','ipnl(2,2)',
                't_intercept','beta2','beta3','beta4']

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
    
def read_visinfo(fname):
    if fname is None:
        return None
        
    colnames = ['bfevis(-2,-2)','bfevis(-1,-2)','bfevis(0,-2)','bfevis(1,-2)','bfevis(2,-2)',
                'bfevis(-2,-1)', 'bfevis(-1,-1)','bfevis(0,-1)','bfevis(1,-1)','bfevis(2,-1)',
                'bfevis(-2,0)','bfevis(-1,0)','bfevis(0,0)','bfevis(1,0)','bfevis(2,0)',
                'bfevis(-2,1)','bfevis(-1,1)','bfevis(0,1)','bfevis(1,1)','bfevis(2,1)',
                'bfevis(-2,2)','bfevis(-1,2)','bfevis(0,2)','bfevis(1,2)','bfevis(2,2)',
                'phi(-2,-2)','phi(-1,-2)','phi(0,-2)','phi(1,-2)','phi(2,-2)',
                'phi(-2,-1)', 'phi(-1,-1)','phi(0,-1)','phi(1,-1)','phi(2,-1)',
                'phi(-2,0)','phi(-1,0)','phi(0,0)','phi(1,0)','phi(2,0)',
                'phi(-2,1)','phi(-1,1)','phi(0,1)','phi(1,1)','phi(2,1)',
                'phi(-2,2)','phi(-1,2)','phi(0,2)','phi(1,2)','phi(2,2)',
                'omega','Cxx','Cxy','Cyy','Ie','p2iter']

    table = pd.read_table(fname,delim_whitespace=True,index_col=(0,1),names=colnames,comment='#')

    # Adjust units as appropriate to this plot
    table['bfevisNN'] = 1e6*(table['bfevis(0,1)']+table['bfevis(0,-1)']
                    + table['bfevis(1,0)']+table['bfevis(-1,0)'])/4.0
    table['phiNN'] = 1e6*(table['phi(0,1)']+table['phi(0,-1)']
                    + table['phi(1,0)']+table['phi(-1,0)'])/4.0
    table['phi(0,0)'] *= 1e6
    table['omega'] *= 100
    return table

font = {'family' : 'normal',
        'weight' : 'regular',
        'size'   : 13}

matplotlib.rc('font', **font)

prefix = '/users/PCON0003/cond0088/Projects/detectors/sw_outputs/PaperIV_chargediffusion/'
summary_files = [prefix+'anna_20663vis_fid1_flatspaperiii_summary.txt',
                 prefix+'ami_modmask_20828vis_fid1_summary.txt',
                 prefix+'anna_20829vis_fid1_summary.txt']

visinfo_files = [prefix+'anna_20663vis_fid1_flatspaperiii_visinfo.txt',
                 prefix+'ami_modmask_20828vis_fid1_visinfo.txt',
                 prefix+'anna_20829vis_fid1_visinfo.txt']

# Set up figure
ylabels = ['SCA 20663',
           'SCA 20828',
           'SCA 20829']
           
divisions = [1,3,5]
axis_names = [r'$\alpha_V$',
              r'$\alpha_H$',
              r'$\beta_2g$',
              r'$\beta_3g^2$',
              r'$\beta_4g^3$',
              r'$g$',
              r'$[K^2a+KK^I]_{0,0,\mathrm{\;vis}}$',
              r'$[K^2a+KK^I]_{<1,0>,\mathrm{\;vis}}$',
              r'$\omega$',
              r'$C_{xx}$',
              r'$C_{yy}$',
              r'$C_{xy}$']
units = ['%',
         '%',
         r'$10^6\times$DN$^{-1}$',
         r'$10^{10}\times$DN$^{-2}$',
         r'$10^{15}\times$DN$^{-3}$',
         r'e/DN',
         r'ppm/e',
         r'ppm/e',
         r'%',
         r'pix$^2$',
         r'pix$^2$',
         r'pix$^2$']

xdata_labels = ['alphaV','alphaH','beta2','beta3','beta4',
                'gain_alphabeta','ipnl(0,0)','ipnlNN','omega','Cxx','Cyy','Cxy']

colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:(divisions[1]-divisions[0])]

fsize = (2*len(axis_names)+5,6) 

fig = plt.figure(figsize=fsize)
grid = plt.GridSpec(16,len(axis_names))
axes=[]

for i in range(len(axis_names)):
    ax = fig.add_subplot(grid[:-2,i])
    axes.append(ax)


summary_tables = [read_summary(f) for f in summary_files]
visinfo_tables = [read_visinfo(f) for f in visinfo_files]

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
        # standard summary values
        if i<len(axes)-4:
            mean = np.mean(summary_tables[j][xdata_labels[i]])
            stdev = np.std(summary_tables[j][xdata_labels[i]])
            error = stdev#/len(summary_tables[j]**0.5)
        
        #visinfo values
        else:
            mean = np.mean(visinfo_tables[j][xdata_labels[i]])
            stdev = np.std(visinfo_tables[j][xdata_labels[i]])
            error = stdev#/len(visinfo_tables[j]**0.5 
        
        #clr = colors[j%(divisions[1]-divisions[0])]
            
        ax.errorbar([mean],[j],xerr=[error],capsize=8.0,markersize=10,marker='o')#,color=clr)
        #if j==0:
        #    for k, div in enumerate(divisions):
        #        ax.axhline(y=div+0.5,color='k',linewidth=0.75,linestyle='-')

        #if j%(divisions[1]-divisions[0])==0:
        #    if j==0 or j==2:
        #        ax.fill_betweenx([j-0.5,j+divisions[1]-divisions[0]-0.5],[mean-error,mean-error],
        #                         [mean+error,mean+error], color = 'grey', alpha = 0.5)
        #    elif j==4:
        #        ax.fill_betweenx([j-0.5,j+divisions[2]-divisions[1]-0.5],[mean-error,mean-error],
        #                         [mean+error,mean+error], color = 'grey', alpha = 0.5)
            

plt.tight_layout()

plt.savefig('paramplot_vis.pdf',format='pdf')
plt.show()

