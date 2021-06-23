import numpy as np
import matplotlib
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')


font = {'family' : 'normal',
        'weight' : 'regular',
        'size'   : 13}

matplotlib.rc('font', **font)
use_cmap = 'gnuplot'

t='output/anna_20829vis_fid1_visinfo.npy'
vis_out_data=np.load(t)


QYomega=vis_out_data[:,:,50]*100
Ie=vis_out_data[:,:,54]/1000

num_bins = 30
F = plt.figure(figsize=(8,6))
S = F.add_subplot(2,2,1)
S.hist(QYomega.ravel(),bins=np.linspace(0, 10, num=num_bins))
S.set_xlabel(r'$\omega$ [%]')

S = F.add_subplot(2,2,2)
S.hist(Ie.ravel(),bins=np.linspace(0.80, 1.10, num=num_bins))
S.set_xlabel(r'$I_e$ [ke/s]')

S = F.add_subplot(2,2,3)
S.hist(vis_out_data[:,:,55].ravel(),bins=np.linspace(0, 100, num=num_bins))
S.set_xlabel(r'Number of iterations')

S = F.add_subplot(2,2,4)
S.hist(vis_out_data[:,:,51].ravel(), num_bins, histtype='step', label=r'$C_{11}$', linewidth=1.5, linestyle='-')
S.hist(vis_out_data[:,:,52].ravel(), bins=np.linspace(-0.1, 0.1, num=num_bins), histtype='step', label=r'$C_{12}$', linewidth=1.5, linestyle='--')
S.hist(vis_out_data[:,:,53].ravel(), num_bins, histtype='step', label=r'$C_{22}$', linewidth=1.5, linestyle='-.')
S.set_xlabel(r'Charge diffusion component [pix$^2$]')
S.legend(loc='upper right', fontsize=12,frameon=False)
F.set_tight_layout(True)
F.savefig('anna_20829vis_fid1_vis_hist.pdf', bbox_inches='tight')

plt.close(F)


dx=128
dy=128

F = plt.figure(figsize=(16,9.5))
S = F.add_subplot(2,3,1)
S.set_title(r'$\omega$ [%]')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,50]*100, cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(2,3,2)
S.set_title(r'$I_e$ [ke/s]')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
#S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,54]/1000, cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(2,3,3)
S.set_title(r'Number of iterations')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
#S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,55], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(2,3,4)
S.set_title(r'$C_{11}$ [pix$^2$]')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,51], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(2,3,5)
S.set_title(r'$C_{12}$ [pix$^2$]')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
#S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,52], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(2,3,6)
S.set_title(r'$C_{22}$ [pix$^2$]')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
#S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,53], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

# F.set_tight_layout(True)
F.savefig('anna_20829vis_fid1_vis_matrices.pdf', bbox_inches='tight')
plt.close(F)

F = plt.figure(figsize=(11, 15))
S = F.add_subplot(3,2,1)
S.set_title(r'$\omega$ [%]')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,50]*100, cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(3,2,2)
S.set_title(r'$I_e$ [ke/s]')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
#S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,54]/1000, cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(3,2,3)
S.set_title(r'Number of iterations')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,55], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(3,2,4)
S.set_title(r'$C_{11}$ [pix$^2$]')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
# S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,51], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(3,2,5)
S.set_title(r'$C_{12}$ [pix$^2$]')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,52], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

S = F.add_subplot(3,2,6)
S.set_title(r'$C_{22}$ [pix$^2$]')
S.set_xlabel('Super pixel X/{:d}'.format(dx))
#S.set_ylabel('Super pixel Y/{:d}'.format(dy))
im = S.imshow(vis_out_data[:,:,53], cmap=use_cmap, origin='lower')
F.colorbar(im, orientation='vertical')

# F.set_tight_layout(True)
F.savefig('anna_20829vis_fid1_vis_matrices_3x2.pdf', bbox_inches='tight')
plt.close(F)
