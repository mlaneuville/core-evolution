#!/usr/local/bin/python
# Time-stamp: <2016-01-28 17:17:08 marine>
# Project : Thermal evolution of stratified core
# Subproject : plot output.
# Author : Matthieu Laneuville, Marine Lasbleis

import matplotlib.pyplot as plt
import numpy as np
import sys


fig_folder = "fig/"
basename = "output-0.txt"

if len(sys.argv)>1:
    N = int(sys.argv[1])
    print N, "profiles to display."
    if len(sys.argv)==3:
        basename = sys.argv[2]+"-"+basename
else: N = 2


# prepare data
# TODO: grid_size should be read from file in future version
data = np.genfromtxt(basename)
grid_size = 200
num_tstep = data.shape[0]/grid_size
print "(grid_size, num_tstep) = (%d,%d)" % (grid_size, num_tstep)
print

# read data from file
radius = data[:,0]
radius = np.reshape(radius, (num_tstep, grid_size))
time = data[:,1]
time = np.reshape(time, (num_tstep, grid_size))
temperature = data[:,2]
temperature = np.reshape(temperature, (num_tstep, grid_size))
conductivity = data[:,3]
conductivity = np.reshape(conductivity, (num_tstep, grid_size))
convect= data[:,4]
convect = np.reshape(convect, (num_tstep, grid_size))

# check where the "convect" field first changes from 0 to 1
# TODO: this does not catch multiple boundaries!
convective_boundary = np.where(np.diff(convect) == 1)[1]*1./grid_size

# FIG1. time, radius temperature map
fig, ax = plt.subplots(figsize=(8,6))
heatmap = ax.imshow(temperature, 
                    extent=[radius[0,0], radius[0,-1], int(time[0,0]), int(time[-1,0])], 
                    aspect='auto', origin='lower')
cbar = fig.colorbar(heatmap)
plt.plot(convective_boundary, time[:,0], 'k', lw=2)
plt.ylim(time[0,0], time[-1,0])
plt.xlabel("Radius [-]")
plt.ylabel("Time [Ma]")
plt.savefig(fig_folder+"2D-temperature-map.eps", format='eps', bbox_inches='tight')
plt.close()

# FIG2. selection of temperature and conductivity profiles
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    
for number in range(N):

    idx = int(number*num_tstep*1./(N-1))
    idx = min(idx, num_tstep-1)
    print "...", idx

    rad, temp, cond = [], [], []

    rad = radius[idx, :]
    temp = temperature[idx, :]
    cond = conductivity[idx, :]

    LABEL = "Time: "+"{:.0f}".format(time[idx,0])+" Ma"
    ax1.plot(rad, temp, label=LABEL)
    ax2.plot(rad, cond, label=LABEL)

ax1.set_ylabel('Temperature')
ax2.set_ylabel('Conductivity')
ax2.set_xlabel('Radius')
ax1.grid()
ax2.grid()
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
ax1.legend(loc=0)
plt.savefig(fig_folder+"1D-profiles.eps", format='eps', bbox_inches='tight')
plt.close()
