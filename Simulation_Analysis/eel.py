#!/usr/bin/env python
# coding: utf-8

# Import necessary libraries.
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

# Define a function which calculates the end to end length of each individual chain in a specified universe.

def eel_his(u, sid):
    
    eelhis = []
    
    for ts in u.trajectory:
        position_1 = u.select_atoms(f"segid {sid}").positions[0]
        position_2 = u.select_atoms(f"segid {sid}").positions[-1]
        eelhis.append(mda.lib.distances.distance_array(position_1, position_2)[0][0])
    return eelhis

# Set simulation-specific variables specifying random seed used and dcd files to be analyzed.
randseed = '11111'
psf_protein_only = 'spider_HA1B1A1B1x2134_silkworm_A1B1x0_protein_only_for_visualization.psf'
equil_dcd = f'equil_{randseed}_unwrap.dcd'
shear_dcd = f'shear_{randseed}_unwrap.dcd'
pull_dcd = f'pull_{randseed}_unwrap.dcd'
relax_dcd = f'equil_after_shear_{randseed}_unwrap.dcd'
tensile_dcd = f'stretch_{randseed}_unwrap.dcd'

# For shear simulation.
# Load all dcd files and protein only psf into an MDAnalysis universe.
u = mda.Universe(f'{psf_protein_only}', 
                [f'{equil_dcd}',f'{shear_dcd}',f'{relax_dcd}',f'{tensile_dcd}'],in_memory=True)

# For pull simulation. Uncomment the below universe command and comment the above universe command.
# Load all dcd files and protein only psf into an MDAnalysis universe.
# u = mda.Universe(f'{psf_protein_only}', 
                # [f'{equil_dcd}',f'{pull_dcd}',f'{relax_dcd}',f'{tensile_dcd}'],in_memory=True)

# Calculate end to end length using the function eel_his.
nofr = len(u.trajectory) # Number of frames in the trajectory
# nochain = int(u.atoms[-1].segid) # Number of chains in the trajectory
nochain = 100 # Reset nochain to only analyze the first 10 chains.
len_his = np.zeros([nochain, nofr]) # History of all chains in the trajectory. Reset to only analyze the first 10 chains.
for ss in range(1,nochain+1):
    len_his[ss-1] = eel_his(u, ss)
    print(ss)
pkl.dump(len_his, open(f'eel_his.pkl','wb'))

# Compute mean and standard deviation of end to end length at each frame in the simulation then plot.
mean = np.mean(len_his, axis=0)

std = np.std(len_his, axis=0)

fig = plt.figure(figsize=[4,3])
ax = fig.add_subplot(111)
plt.errorbar(np.arange(nofr), mean, std, linestyle='None', marker='^')
plt.xlabel('frame no.')
plt.ylabel(f"End to End Length")
plt.tight_layout()
plt.savefig(f'eel_his.png', dpi=200)
plt.close()




