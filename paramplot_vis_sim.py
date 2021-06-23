import numpy as np
import matplotlib
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.ticker import FormatStrFormatter

sim_gain_od = 1.5
sim_gain = 1.73
sim_alpha = 1.38
sim_alphaH = 1.38
sim_alphaV = 1.38
noipc_sim_alphaH = 0.0
noipc_sim_alphaV = 0.0
#sim_beta2 = -1.5725*1e-6
#sim_beta3 = 1.9307e-5*1e-6
#sim_beta4 = -1.4099e-10*1e-6
sim_beta2 = 0.0
sim_beta3 = 0.0
sim_beta4 = 0.0
sim_ipnl_od = -1.372
sim_ipnlNN_od = 0.2840
sim_ipnl_ir = -1.7779
sim_ipnlNN_ir = 0.3001
sim_ipnl_vis = -2.6777
sim_ipnlNN_vis = 0.3033
noipc_sim_ipnl_ir = -2.0356
noipc_sim_ipnlNN_ir = 0.3886
noipc_sim_ipnl_vis = -3.0462
noipc_sim_ipnlNN_vis = 0.4176
#sim_phi = 0
#sim_phiNN = 0
sim_omega = 0.08*100 # To units of %
sim_Cxx = 0.04
sim_Cxy = 0.02
sim_Cyy = 0.04
dummy = 100.0

sim_truevals = [sim_alphaV,sim_alphaH,
                #sim_beta2*sim_gain*-1e6,sim_beta3*sim_gain**2*-1e10,sim_beta4*sim_gain**3*-1e15,
                sim_gain,sim_ipnl_ir,sim_ipnlNN_ir,
                sim_ipnl_vis,sim_ipnlNN_vis,sim_omega,
                sim_Cxx,sim_Cyy,sim_Cxy]
sim_truevals_noipc = [noipc_sim_alphaV, noipc_sim_alphaH,
                      dummy, noipc_sim_ipnl_ir, noipc_sim_ipnlNN_ir,
                     noipc_sim_ipnl_vis, noipc_sim_ipnlNN_vis,
                     dummy, dummy, dummy, dummy]
sim_truevals_od = [dummy, dummy,
                  sim_gain_od,sim_ipnl_od,sim_ipnlNN_od,
                  sim_ipnl_od,sim_ipnlNN_od,dummy, dummy, dummy, dummy]
sim_ylabels = ['simulations, non-zero cxy','simulations, diff ir/vis bfe']
 
 
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
    table['ipnlNN_ir'] = 1e6*(table['ipnl(0,1)']+table['ipnl(0,-1)']
                    + table['ipnl(1,0)']+table['ipnl(-1,0)'])/4.0
    table['ipnl(0,0)_ir'] = table['ipnl(0,0)']*1e6
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
    table['bfevis(0,0)'] *= 1e6
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

prefix = '/users/PCON0003/cond0088/Projects/detectors/sw_outputs/PaperIV_chargediffusion/sim_out/'
#summary_files = [prefix+'vis-sim.offdiagp_feb3_summary.txt',
summary_files = [prefix+'vis-sim.offdiagp_jun8_summary.txt',
	             prefix+'vis-sim.diffbfe_mar30_summary.txt',
                 prefix+'vis-sim.diffbfe_16x16_summary.txt',
	             prefix+'vis-sim.diffbfe_apr7withlinipc_summary.txt']

#visinfo_files = [prefix+'vis-sim.offdiagp_feb3_visinfo.txt',
visinfo_files = [prefix+'vis-sim.offdiagp_jun8_visinfo.txt',
	             prefix+'vis-sim.diffbfe_mar30_visinfo.txt',
                 prefix+'vis-sim.diffbfe_16x16_visinfo.txt',
	             prefix+'vis-sim.diffbfe_apr7withlinipc_visinfo.txt']

# Set up figure
ylabels = ['Sim, \nnon-zero Cxy',
           'Sim, \ndiff ir/vis BFE, \nno IPC',
           'Sim, \ndiff ir/vis BFE, \nno IPC 16x16',
           'Sim, \ndiff ir/vis BFE, \n with IPC']
           
divisions = [1,3,5] # Where the horizontal lines get drawn
axis_names = [r'$\alpha_V$',
              r'$\alpha_H$',
#              r'$\beta_2g$',
#              r'$\beta_3g^2$',
#              r'$\beta_4g^3$',
              r'$g$',
              r'$[K^2a]_{0,0,ir}$',
              r'$[K^2a]_{<1,0>,ir}$',
              r'$[K^2a]_{0,0,vis}$',
              r'$[K^2a]_{<1,0>,vis}$',
              r'$\omega$',
              r'$C_{xx}$',
              r'$C_{yy}$',
              r'$C_{xy}$']
units = ['%',
         '%',
#         r'$10^6\times$DN$^{-1}$',
#         r'$10^{10}\times$DN$^{-2}$',
#         r'$10^{15}\times$DN$^{-3}$',
         r'e/DN',
         r'ppm/e',
         r'ppm/e',
         r'ppm/e',
         r'ppm/e',
         r'%',
         r'pix$^{2}$',
         r'pix$^{2}$',
         r'pix$^{2}$']

xdata_labels = ['alphaV','alphaH',#'beta2','beta3','beta4',
                'gain_alphabeta','ipnl(0,0)_ir','ipnlNN_ir','bfevis(0,0)','bfevisNN','omega','Cxx','Cyy', 'Cxy']

colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:(divisions[1]-divisions[0])]
colors.append('Gray')
fsize = (2*len(axis_names),6) 

fig = plt.figure(figsize=fsize)
grid = plt.GridSpec(16,len(axis_names))
axes=[]
#sim_axes = []
for i in range(len(axis_names)):
    ax = fig.add_subplot(grid[:-2,i])
    axes.append(ax)
#
#    sax = fig.add_subplot(grid[-2:,i])
#    sim_axes.append(sax)

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
    #ax.xaxis.set_major_formatter('%3.1f')  # This not working yet
    ax.xaxis.labelpad=15

    for j,run in enumerate(ylabels):
        # standard summary values
        if i<len(axes)-6: # this changes depending on which visible params we want to show
            mean = np.mean(summary_tables[j][xdata_labels[i]])
            stdev = np.std(summary_tables[j][xdata_labels[i]])
            error = stdev#/len(summary_tables[j]**0.5)
        
        #visinfo values
        else:
            mean = np.mean(visinfo_tables[j][xdata_labels[i]])
            stdev = np.std(visinfo_tables[j][xdata_labels[i]])
            error = stdev#/len(visinfo_tables[j]**0.5 
        
        clr = colors[j%(divisions[1]-divisions[0])]
        # This makes the script way messier than Jenna's original code
        if (j==0):
            ax.errorbar([mean],[j],xerr=[error],capsize=8.0,markersize=10,marker='o',color=clr)
        elif ((j>0)&(j<3)):
            ax.errorbar([mean],[j],xerr=[error],capsize=8.0,markersize=10,marker='o',color=colors[2])
        else:
            ax.errorbar([mean],[j],xerr=[error],capsize=8.0,markersize=10,marker='o',color=clr)
        #if j==0:
        #    for k, div in enumerate(divisions):
        #        ax.axhline(y=div+0.5,color='k',linewidth=0.75,linestyle='-')

        if j%(divisions[1]-divisions[0])==0:
            ax.axvline(x=sim_truevals[i],linestyle='--',color='k')
            if ((i==2)|(i>6)):
                pass
            else:
                ax.axvline(x=sim_truevals_noipc[i],linestyle=':',color='Gray')
            #if j==0 or j==2:
            #    ax.fill_betweenx([j-0.5,j+divisions[1]-divisions[0]-0.5],[mean-error,mean-error],
            #                     [mean+error,mean+error], color = 'grey', alpha = 0.5)
            #elif j==4:
            #    ax.fill_betweenx([j-0.5,j+divisions[2]-divisions[1]-0.5],[mean-error,mean-error],
            #                     [mean+error,mean+error], color = 'grey', alpha = 0.5)
            
sim_tables = summary_tables[-1*len(sim_ylabels):]

#for i,sax in enumerate(sim_axes):
#    sax.set_ylim(len(sim_ylabels)-0.5,-0.5)
#    sax.ticklabel_format(axis='x', style='sci', scilimits=(-3,3))
#    sax.set_yticklabels([])
#    
#    if i==0:
#        sax.set_yticks([n for n in range(len(sim_ylabels))])
#        sax.set_yticklabels(sim_ylabels)
#        sax.set_ylim(len(sim_ylabels)-0.5,-0.5)
#
#    sax.tick_params(axis='y',which='both',left=False)
#    
#    for j,run in enumerate(sim_ylabels):
#	mean = np.mean(sim_tables[j][xdata_labels[i]])
#	stdev = np.std(sim_tables[j][xdata_labels[i]]) # Possible error in paramplot.py?????
#       error = stdev#/len(tables[j]**0.5)
#	sax.errorbar([mean],[j],xerr=[error],capsize=8.0,markersize=10,marker='o',
#                    color=colors[j%(divisions[1]-divisions[0])])
#
#        if j%(divisions[1]-divisions[0])==0:
#	    sax.axvline(x=sim_truevals[i],linestyle='--',color='k')
#            #sax.fill_betweenx([j-0.5,j+len(sim_ylabels)-0.5],[mean-stdev,mean-stdev],
#            #                 [mean+stdev,mean+stdev], color = 'grey', alpha = 0.5)

plt.tight_layout()

plt.savefig('paramplot_vis_sim.pdf',format='pdf')
plt.show()

