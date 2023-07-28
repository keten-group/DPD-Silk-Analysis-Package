import numpy as np

import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('TkAgg')
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as plticker
plt.rcParams["font.family"] = "Arial"
from matplotlib import cm

import os
####-----------------read data
# time in fs
sys_names = ['HA1B1A1B1', 'HA1B3', 'HA3B1']

num_beta={}
ave_crystal_size={}
num_conn={}
conduct={}

for ss in range(len(sys_names)):

	with open(f'{sys_names[ss]}_simulation/equil.txt') as f:
	    lines = f.readlines()
	data = np.zeros([len(lines),2])
	ii = 0
	for line in lines:
		data[ii] = line.strip('\n').split(' ')
		ii += 1
	totframe = int(len(data)//3)
	num_beta_equil = data[:totframe]
	ave_crystal_size_equil = data[totframe:2*totframe]
	num_conn_equil = data[2*totframe:]

	with open(f'{sys_names[ss]}_simulation/shear.txt') as f:
	    lines = f.readlines()
	data = np.zeros([len(lines),2])
	ii = 0
	for line in lines:
		data[ii] = line.strip('\n').split(' ')
		ii += 1
	totframe = int(len(data)//3)
	num_beta_shear = data[:totframe]
	ave_crystal_size_shear = data[totframe:2*totframe]
	num_conn_shear = data[2*totframe:]

	num_beta_shear[:,0] += num_beta_equil[:,0][-1]
	ave_crystal_size_shear[:,0] += ave_crystal_size_equil[:,0][-1]
	num_conn_shear[:,0] += num_conn_equil[:,0][-1]

	num_beta[sys_names[ss]]=np.vstack([num_beta_equil,num_beta_shear])
	ave_crystal_size[sys_names[ss]]=np.vstack([ave_crystal_size_equil,ave_crystal_size_shear])
	num_conn[sys_names[ss]]=np.vstack([num_conn_equil,num_conn_shear])

	# conductance
	data = []
	for ie in range(1,8):
		thisdir = f'{sys_names[ss]}_simulation/equil_evolve_{ie}/conductance.txt'
		if os.path.exists(thisdir):
			with open(thisdir) as f:
				    lines = f.readlines()
			data.append(lines[0])
		else:
			data.append(0)
	for ishear in range(15):
		thisdir = f'{sys_names[ss]}_simulation/shear_evolve_{ishear}/conductance.txt'
		if os.path.exists(thisdir):
			with open(thisdir) as f:
				    lines = f.readlines()
			data.append(lines[0])
		else:
			data.append(0)
	data = np.array(data).astype(float)
	conduct[sys_names[ss]]=np.vstack([num_conn[sys_names[ss]][:,0],data]).T

####---------------plot

##### Number of Beta Crystals
fig = plt.figure(figsize=[6, 4])
ax = fig.add_subplot(111)

for name in sys_names:
	plt.plot(num_beta[name][:,0], num_beta[name][:,1], '.-', markersize=10, label=name)

ymin, ymax = ax.get_ylim()
plt.vlines(num_beta_equil[:,0][-1], ymin, ymax, linestyle='--', linewidth=1.0, color='k')
plt.legend(loc='best')
plt.xlim(0, )
plt.ylim(0, )
plt.xlabel('Simulation Time [fs]', fontsize=11)
plt.ylabel('Number of Beta Crystals', fontsize=11)
plt.tight_layout()
fig.savefig(f'num_beta.png', format="png", dpi=200)
plt.close()


##### Average Crystal Size
fig = plt.figure(figsize=[6, 4])
ax = fig.add_subplot(111)

for name in sys_names:
	plt.plot(ave_crystal_size[name][:,0], ave_crystal_size[name][:,1], '.-', markersize=10, label=name)

ymin, ymax = ax.get_ylim()
plt.vlines(ave_crystal_size_equil[:,0][-1], ymin, ymax, linestyle='--', linewidth=1.0, color='k')
plt.legend(loc='best')
plt.xlim(0, )
plt.ylim(0, )
plt.xlabel('Simulation Time [fs]', fontsize=11)
plt.ylabel('Average Crystal Size', fontsize=11)
plt.tight_layout()
fig.savefig(f'ave_crystal_size.png', format="png", dpi=200)
plt.close()

##### Number of Connections
fig = plt.figure(figsize=[6, 4])
ax = fig.add_subplot(111)

for name in sys_names:
	plt.plot(num_conn[name][:,0], num_conn[name][:,1], '.-', markersize=10, label=name)

ymin, ymax = ax.get_ylim()
plt.vlines(num_conn_equil[:,0][-1], ymin, ymax, linestyle='--', linewidth=1.0, color='k')
plt.legend(loc='best')
plt.xlim(0, )
plt.ylim(0, )
plt.xlabel('Simulation Time [fs]', fontsize=11)
plt.ylabel('Number of Beta Connections', fontsize=11)
plt.tight_layout()
fig.savefig(f'num_conn.png', format="png", dpi=200)
plt.close()

##### Network conductance
fig = plt.figure(figsize=[6, 4])
ax = fig.add_subplot(111)

for name in sys_names:
	plt.plot(conduct[name][:,0], conduct[name][:,1], '.-', markersize=10, label=name)

ymin, ymax = ax.get_ylim()
plt.vlines(num_conn_equil[:,0][-1], ymin, ymax, linestyle='--', linewidth=1.0, color='k')
plt.legend(loc='best')
plt.xlim(0, )
plt.ylim(0, )
plt.xlabel('Simulation Time [fs]', fontsize=11)
plt.ylabel('Conductance', fontsize=11)
plt.tight_layout()
fig.savefig(f'conductance.png', format="png", dpi=200)
plt.close()
