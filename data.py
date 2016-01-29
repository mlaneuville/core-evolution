#!/usr/local/bin/python
# Time-stamp: <2016-01-29 12:09:28 marine>
# Project : Thermal evolution of stratified core
# Subproject : Read (and correct/interpolate/etc.) data sets given by G.H. 
# Author : Marine Lasbleis


import matplotlib.pyplot as plt
import numpy as np
import sys



def homogeneous_sample_data(basename="earth", Npoints=200):


    filename = basename+".dat"
    output_file = basename+"_resampled.txt"

    fig, ax = plt.subplots(3, 2, sharex=True)
    
    fig.suptitle("Data and fit", fontsize=14)
    ax[0,0].set_ylabel('Temperature (K)')
    ax[1,0].set_ylabel('Pressure (Pa)')
    ax[2,0].set_ylabel('Gravity (m/s**2)')
    ax[2,0].set_xlabel('Radius (km)')
    ax[0,1].set_ylabel('(Temperature -T0)/Radius**2 (K/m**2)')
    ax[1,1].set_ylabel('(Pressure -P0)/Radius**2) (Pa/m**2)')
    ax[2,1].set_ylabel('Gravity/Radius (1./m/s**2)')
    ax[2,1].set_xlabel('Radius (km)')
    
    data = np.genfromtxt(filename)
    print "====="
    print "Number of points, number of columns : ", np.shape(data)
    print "Data is R (km) T (K) g (m/s**2) P (Pa)"
    print "====="
    Radius = data[:,0]*1e3 #in m
    Temperature = data[:,1] #in K
    Gravity = data[:,2] #in m/s**2
    Pressure = data[:,3] #in Pa

    if abs((Temperature[0]-Temperature[-1])/Temperature[-1]) < 0.01:
        print "Please check the data. Your temperature profile is almost constant."
        print "The data file is ", name
        print "Temperature profile goes from ", Temperature[0], "K to ", Temperature[-1], "K."
        print "====="
        

    ax[0,0].plot(Radius/1e3, Temperature, '+r')
    ax[1,0].plot(Radius/1e3, Pressure, '+r')
    ax[2,0].plot(Radius/1e3, Gravity, '+r')

    # functions Temperature and Pressure are fitted with polynomes
    # Constraints are : value at 0 (T[0] and P[0]) and tangent at 0 (=0)
    # Gravity is approximated by g=g0 r (linear in radius)
    A0 = Temperature[0]
    B0 = Pressure[0]
    A1 = 0. #tangent at x=0

    
    ax[0,1].plot(Radius[1:]/1e3, (Temperature[1:] -(A0))/Radius[1:]**2., 'r+')
    ax[1,1].plot(Radius[1:]/1e3, (Pressure[1:] -(B0))/Radius[1:]**2., 'r+')
    ax[2,1].plot(Radius[1:]/1e3, Gravity[1:]/Radius[1:], 'r+')

    polyn_temp = np.polyfit(Radius[1:], (Temperature[1:] -(A0))/Radius[1:]**2., 4)
    polyn_press = np.polyfit(Radius[1:], (Pressure[1:] -(B0))/Radius[1:]**2., 4)
    polyn_grav = np.polyfit(Radius[1:], Gravity[1:]/Radius[1:], 0)

    #new radius repartition (resampling with constant dr)
    r = np.linspace(Radius[0], Radius[-1], Npoints)
    temp = np.polyval(polyn_temp, r)*r**2+A0
    press = np.polyval(polyn_press, r)*r**2+B0
    grav = np.polyval(polyn_grav, r)*r
    
    ax[0,0].plot(r/1e3, temp)
    ax[1,0].plot(r/1e3, press)
    ax[2,0].plot(r/1e3, grav)
    ax[0,1].plot(r/1e3, np.polyval(polyn_temp, r))
    ax[1,1].plot(r/1e3, np.polyval(polyn_press, r))
    ax[2,1].plot(r/1e3, np.polyval(polyn_grav, r))
    ax[2,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    

    data_resample = np.concatenate((np.array([r]).T, np.array([temp]).T, np.array([press]).T, np.array([grav]).T), axis=1)
    np.savetxt(output_file, data_resample, header="radius (m), temperature (K)")

    print "red '+': data points, blue line: fit and resample"
    plt.show()
    
    return



if __name__ == '__main__':


    # name of the data file (default: earth.dat)
    if len(sys.argv)>1:
        name = sys.argv[1]
        print "data file: ", name
    else: name = "earth"

    # number of points in the new profile
    if len(sys.argv)==3:
        Npoints = int(sys.argv[2])
    else: Npoints = 200


    homogeneous_sample_data(name, Npoints)
