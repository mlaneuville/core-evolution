#!/usr/local/bin/python
# Time-stamp: <2016-02-16 09:54:51 marine>
# Project : Thermal evolution of stratified core
# Subproject : plot output.
# Author : Matthieu Laneuville, Marine Lasbleis

import matplotlib.pyplot as plt
import numpy as np
import sys



def read_data_from_file(basename, grid_size=200, out_folder="out/", constantradius="True"):

    # prepare data
    # TODO: grid_size should be read from file in future version
    data = np.genfromtxt(out_folder+basename+"-output.txt")
    num_tstep = data.shape[0]/grid_size
    print "(grid_size, num_tstep) = (%d,%d)" % (grid_size, num_tstep)
    print

    # read data from file
    # TODO: loading everything like this sounds like a terrible idea; can we improve? I don't think so...
    # But I changed the organisation and add a function, so that the "data" is actually only a local variable.
    # TODO : maybe cut here the data if only the profiles are needed? (only Fig2)
    
    radius = data[:,0]
    radius = np.reshape(radius, (num_tstep, grid_size))
    radius = np.array(radius[0,:])
    time = data[:,1]
    time = np.reshape(time, (num_tstep, grid_size))
    temperature = data[:,2]
    temperature = np.reshape(temperature, (num_tstep, grid_size))
    conductivity = data[:,3]
    conductivity = np.reshape(conductivity, (num_tstep, grid_size))
    adiabat = data[:,4]
    adiabat = np.reshape(adiabat, (num_tstep, grid_size))
    qcmb = data[:,5]
    qcmb = np.reshape(qcmb, (num_tstep, grid_size))
    convect = data[:,6]
    convect = np.reshape(convect, (num_tstep, grid_size))


    return num_tstep, radius, time, temperature, conductivity, adiabat, qcmb, convect


if __name__ == '__main__':

    fig_folder = "fig/"
    out_folder = "out/"

    if len(sys.argv) != 3:
        print "You have to specify the number of profiles and the name of the run."
        print "Usage: python "+sys.argv[0]+" <run_name> <num_profiles>"
        sys.exit()
    else:
        basename = sys.argv[1]
        N = int(sys.argv[2])
        print basename, ", ", N, "profiles to display."

   
    grid_size = 200
    num_tstep, radius, time, temperature, conductivity, adiabat, qcmb, convect = read_data_from_file(basename, grid_size, out_folder)


    # check where the "convect" field first changes from 0 to 1
    # TODO: this does not catch multiple boundaries!
    convective_boundary = radius[(np.where(np.diff(convect) == 1)[1])]#*1.*radius[0,-1]/grid_size


    # FIG1. time, radius temperature map
    fig, ax = plt.subplots(figsize=(12,8))
    heatmap = ax.imshow(temperature, 
                        extent=[radius[0], radius[-1], int(time[0,0]), int(time[-1,0])], 
                        aspect='auto', origin='lower')
    cbar = fig.colorbar(heatmap)
    plt.plot(convective_boundary, time[:,0], 'k', lw=2)
    plt.xlim(radius[0], radius[-1])
    plt.ylim(time[0,0], time[-1,0])
    plt.xlabel("Radius [km]")
    plt.ylabel("Time [Ma]")
    plt.savefig(fig_folder+basename+"-2D-temperature-map.eps", format='eps', bbox_inches='tight')
    plt.close()

    # FIG2. selection of temperature and conductivity profiles
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(10,10))

    for number in range(N):

        idx = int(number*num_tstep*1./(N-1))
        idx = min(idx, num_tstep-1)
        print "...", idx

        rad, temp, cond, ad = [], [], [], []

        rad = radius[:]
        temp = temperature[idx, :]
        cond = conductivity[idx, :]
        ad = adiabat[idx, :]*1000
        # TODO: generalize this
        dtemp = np.append(np.array([0]), np.diff(temp))/(rad[-1]/grid_size)

        LABEL = "Time: "+"{:.0f}".format(time[idx,0])+" Ma"

        ax1.plot(rad, temp, label=LABEL, lw=2)
        ax2.plot(rad, cond/1e-5, label=LABEL, lw=2)

        ax3.plot(rad, -ad, label=LABEL, lw=2)
        ax4.plot(rad, dtemp, label=LABEL, lw=2)
        ax6.plot(rad, -ad, label=LABEL, lw=2, ls='--')
        ax6.plot(rad, dtemp, label=LABEL, lw=2)

    ax5.plot(time[:,0], qcmb[:,0], lw=2)
    ax5.set_xlim(time[0,0], time[-1,0])

    ax1.set_xlim(rad[0], rad[-1])
    ax2.set_xlim(rad[0], rad[-1])
    ax3.set_xlim(rad[0], rad[-1])
    ax4.set_xlim(rad[0], rad[-1])
    ax6.set_xlim(rad[0], rad[-1])

    ax1.set_ylabel('Temperature [K]')
    ax2.set_ylabel('Diffusivity [1e-5 m$^2$/s]')
    ax3.set_ylabel('Adiabatic gradient [K/km]')
    ax4.set_ylabel('Temperature gradient [K/km]')
    ax5.set_ylabel('CMB heat flow [TW]')
    ax6.set_ylabel('Temperature gradient [K/km]')

    ax1.set_xlabel('Radius [km]')
    ax2.set_xlabel('Radius [km]')
    ax3.set_xlabel('Radius [km]')
    ax4.set_xlabel('Radius [km]')
    ax5.set_xlabel('Time [Ma]')
    ax6.set_xlabel('Radius [km]')

    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()
    ax5.grid()
    ax6.grid()


    #ax1.legend(loc=0)
    plt.tight_layout()
    #plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    plt.savefig(fig_folder+basename+"-1D-profiles.eps", format='eps', bbox_inches='tight')
    plt.close()
