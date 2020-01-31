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
            'ipnl(2,2)','t_intercept','beta2g','beta3g2','beta4g3']

    table = pd.read_table(fname,delim_whitespace=True,index_col=(0,1),names=colnames,comment='#')

    # Adjust units 
    for n in colnames:
        if n[:4]=='ipnl':
            table[n] *= 1e6
            
        if n[:5]=='alpha':
            table[n] *= 100.
            
        if n=='q_per_t':
            table[n] *= 0.001
            
        if n=='beta2g':
            table[n] *= -1e6
            
        if n=='beta3g2':
            table[n] *= -1e9
        
        if n=='beta4g3':
            table[n] *= -1e15
            
    # Average quantities

    table['alpha'] = (2*table['alphaH']+2*table['alphaV']+4*table['alphaD'])/8.
    
    table['ipnl<1,0>'] = (table['ipnl(0,1)']+table['ipnl(0,-1)']
                         +table['ipnl(1,0)']+table['ipnl(-1,0)'])/4.0
    table['ipnl<1,1>'] = (table['ipnl(1,1)']+table['ipnl(1,-1)']
                         +table['ipnl(-1,1)']+table['ipnl(-1,-1)'])/4.0
    table['ipnl<2,0>'] = (table['ipnl(0,2)']+table['ipnl(0,-2)']
                         +table['ipnl(2,0)']+table['ipnl(-2,0)'])/4.0
    table['ipnl<2,1>'] = (table['ipnl(2,1)']+table['ipnl(2,-1)']
                         +table['ipnl(-2,1)']+table['ipnl(-2,-1)']
                         +table['ipnl(1,2)']+table['ipnl(1,-2)']
                         +table['ipnl(-1,2)']+table['ipnl(-1,-2)'])/8.0
    table['ipnl<2,2>'] = (table['ipnl(2,2)']+table['ipnl(2,-2)']
                         +table['ipnl(-2,2)']+table['ipnl(-2,-2)'])/4.0
                             
    table['ipnlH'] = (table['ipnl(1,0)']+table['ipnl(-1,0)'])/2.0
    table['ipnlV'] = (table['ipnl(0,1)']+table['ipnl(0,-1)'])/2.0
                             
    
    return table

prefix = '/users/PCON0003/cond0088/Projects/detectors/sw_outputs/'
files_20663 = [prefix+'chris_20663st_summary.txt',
	           prefix+'chris_20663st_128x16_summary.txt',
	           prefix+'chris_20663st-cub_summary.txt',
	           prefix+'chris_20663st-lo_summary.txt',
	           prefix+'chris_20663st-short_summary.txt',
                   prefix+'chris_20663st-med_summary.txt',        
	           prefix+'chris_20663st-e2_summary.txt']

files_20828 = [prefix+'chris_20828st_summary.txt',
	           prefix+'chris_20828st_128x16_summary.txt',
	           prefix+'chris_20828st-cub_summary.txt',
 	           prefix+'chris_20828st-lo_summary.txt',
	           prefix+'chris_20828st-short_summary.txt',
               	   prefix+'chris_20828st-med_summary.txt',
		   prefix+'chris_20828st-e2_summary.txt']
        
files_20829 = [prefix+'chris_20829st_summary.txt',
	           prefix+'chris_20829st_128x16_summary.txt',
	           prefix+'chris_20829st-cub_summary.txt',
	           prefix+'chris_20829st-lo_summary.txt',
	           prefix+'chris_20829st-short_summary.txt',
	           prefix+'chris_20829st-med_summary.txt',
		   prefix+'chris_20829st-e2_summary.txt']


rows = ['q_per_t','gain_alphabeta','alpha','alphaH','alphaV','alphaD',
        'beta2g','beta3g2','beta4g3','ipnl(0,0)','ipnl<1,0>','ipnl<1,1>',
        'ipnl<2,0>','ipnl<2,1>','ipnl<2,2>','ipnlH','ipnlV']
# Print in tex-friendly format
for block in [files_20663,files_20828,files_20829]:
    fid_all = read_summary(block[0])
    fid = np.mean(fid_all)
    fid_sd = np.std(fid_all)

    s128x16_all = read_summary(block[1])
    s128x16 = np.mean(s128x16_all)
    s128x16_sd = np.std(s128x16_all)

    cubic_all = read_summary(block[2])
    cubic = np.mean(cubic_all)
    cubic_sd = np.std(cubic_all)

    lo_all = read_summary(block[3])
    lo = np.mean(lo_all)
    lo_sd = np.std(lo_all)

    short_all = read_summary(block[4])
    short = np.mean(short_all)
    short_sd = np.std(short_all)

    med_all = read_summary(block[5])
    med = np.mean(med_all)
    med_sd = np.std(med_all)

    eps_all = read_summary(block[6])
    eps = np.mean(eps_all)
    eps_sd = np.std(eps_all)
    
    print('Means')
    for r in rows:
	print '&{:0.4f} &{:0.4f} &{:0.4f} &{:0.4f} &{:0.4f} &{:0.4f} &{:0.4f} \\\\'.format(
                                        fid[r],s128x16[r],cubic[r],lo[r],short[r],med[r],eps[r])	
	if r=='beta4g3':
	    print('\n')

    print('\nStandard deviations')
    for r in rows:
        print '&{:0.4f} &{:0.4f} &{:0.4f} &{:0.4f} &{:0.4f} &{:0.4f} &{:0.4f} \\\\'.format(
                        fid_sd[r],s128x16_sd[r],cubic_sd[r],lo_sd[r],short_sd[r],med_sd[r],eps_sd[r])            
    	if r=='beta4g3':
            print('\n')

    print('\n')

