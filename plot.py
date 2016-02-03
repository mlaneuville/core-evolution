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
suffix = ""

if len(sys.argv)>1:
    N = int(sys.argv[1])
    print N, "profiles to display."
    if len(sys.argv)==3:
        basename = sys.argv[2]+"-"+basename
        suffix = "-"+sys.argv[2]
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
adiabat = data[:,4]
adiabat = np.reshape(adiabat, (num_tstep, grid_size))
convect = data[:,5]
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
plt.savefig(fig_folder+"2D-temperature-map"+suffix+".eps", format='eps', bbox_inches='tight')
plt.close()

# FIG2. selection of temperature and conductivity profiles
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    
for number in range(N):

    idx = int(number*num_tstep*1./(N-1))
    idx = min(idx, num_tstep-1)
    print "...", idx

    rad, temp, cond, ad = [], [], [], []

    rad = radius[idx, :]
    temp = temperature[idx, :]
    cond = conductivity[idx, :]
    ad = adiabat[idx, :]*1000
    dtemp = np.append(np.array([0]), np.diff(temp))/(3480e3/grid_size)*1000

    LABEL = "Time: "+"{:.0f}".format(time[idx,0])+" Ma"

    ax1.plot(rad, temp, label=LABEL)
    ax2.plot(rad, cond, label=LABEL)

    ax3.plot(rad, -ad, label=LABEL)
    ax4.plot(rad, dtemp, label=LABEL)

ax1.set_ylabel('Temperature [K]')
ax2.set_ylabel('Conductivity [-]')
ax3.set_ylabel('Adiabatic gradient [K/km]')
ax4.set_ylabel('Temperature gradient [K/km]')

ax1.set_xlabel('Radius [-]')
ax2.set_xlabel('Radius [-]')
ax3.set_xlabel('Radius [-]')
ax4.set_xlabel('Radius [-]')

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()

#ax1.legend(loc=0)
plt.tight_layout()
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
plt.savefig(fig_folder+"1D-profiles"+suffix+".eps", format='eps', bbox_inches='tight')
plt.close()
