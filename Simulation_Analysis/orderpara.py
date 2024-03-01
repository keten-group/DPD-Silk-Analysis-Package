#!/usr/bin/env python
# coding: utf-8

import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np
import pickle as pkl
import os.path
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as plticker

def order_his(u, sid):
    
    shis = []

    for ts in u.trajectory:
        temppos = u.select_atoms(f"segid {sid}").positions[1:]-u.select_atoms(f"segid {sid}").positions[:-1]
        temppos = temppos/np.linalg.norm(temppos,axis=1).reshape([len(temppos),1])
        shis.append(0.5*(3*np.mean((np.inner(temppos, np.array([1,0,0])))**2)-1))

    return np.array(shis)

randseed = '11111'
psf_protein_only = 'spider_HA1B1A1B1x2134_silkworm_A1B1x0_protein_only.psf'
equil_dcd = f'equil_{randseed}_unwrap.dcd' 
pull_dcd = f'pull_{randseed}_unwrap.dcd'


u = mda.Universe(f'{psf_protein_only}', 
                [f'{equil_dcd}',f'{pull_dcd}'],in_memory=True)

nofr = len(u.trajectory) # Number of frames in the trajectory
# nochain = int(u.atoms[-1].segid) # Number of chains in the trajectory
nochain = 100 # Reset nochain to only analyze the first 100 chains.
len_his = np.zeros([nochain, nofr]) # History of all chains in the trajectory. Reset to only analyze the first 100 chains.
for ss in range(1,nochain+1):
    len_his[ss-1] = order_his(u, ss)
    print(ss)

mean = np.mean(len_his, axis=0)

std = np.std(len_his, axis=0)

fig = plt.figure(figsize=[4,3])
ax = fig.add_subplot(111)
plt.errorbar(np.arange(nofr), mean, std, linestyle='None', marker='^')
plt.xlabel('frame no.')
plt.ylabel(f"Herman's Order")
# plt.ylim([-0.005,0.005])
plt.tight_layout()
plt.savefig(f'order_his.png', dpi=200)
plt.close()

pkl.dump(len_his, open(f'order_his.pkl','wb'))



