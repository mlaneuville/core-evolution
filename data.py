#!/usr/local/bin/python
# Time-stamp: <2016-02-01 11:46:13 marine>
# Project : Thermal evolution of stratified core
# Subproject : Read (and correct/interpolate/etc.) data sets given by G.H. 
# Author : Marine Lasbleis


import matplotlib.pyplot as plt
import numpy as np
import sys


### `python data.py` will run over the 4 planets (Earth, Venus, Mars, Vesta)
### `python data.py planet` will run over the planet `planet` only


def homogeneous_sample_data(basename="earth", folder="./dat/", Npoly=1, resample=0, Npoints=20):
    """ Open file folder/basename.dat and transform the data to be used by main.cpp

    parameters:
    - basename, folder: data file (folder/basename.dat will be used. Default is ./dat/earth.dat)
    data file has 4 columns: R (km) T (K) g (m/s**2) P (Pa). 
    - Npoly: (Npoly+1) is the number of free coefficients in the fit (default is 1, which means highest order polynome is 3)
    - resample: decide if a resample file is output or not. (0/False by default)
    - Npoints: number of points for the resampling. Default is 20.
    return: nothing.
    print 1 or 2 file:
    - polynomes_basename: polynomial coefficients
    and if resample=1, another file is also printed:
    - resampled_basename: similar data set than the input ones, but resampled with constant steps in radius.
    (output is SI, so radius in m and pressure in Pa!)
    """

    filename = folder+basename+".dat"
    output_file1 = folder+"polynomes_"+basename+".dat"
    output_file2 = folder+"resampled_"+basename+".dat"
    
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
    print "data file: ", filename
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
    T0 = Temperature[0]
    P0 = Pressure[0]
    A1 = 0. #tangent at x=0

    ax[0,1].plot(Radius[1:]/1e3, (Temperature[1:] -(T0))/Radius[1:]**2., 'r+')
    ax[1,1].plot(Radius[1:]/1e3, (Pressure[1:] -(P0))/Radius[1:]**2., 'r+')
    ax[2,1].plot(Radius[1:]/1e3, Gravity[1:]/Radius[1:], 'r+')
    
    polyn_temp = np.polyfit(Radius[1:], (Temperature[1:] -(T0))/Radius[1:]**2., Npoly)
    polyn_press = np.polyfit(Radius[1:], (Pressure[1:] -(P0))/Radius[1:]**2., Npoly)
    polynomes_temp = np.append(polyn_temp, np.array([0, T0]))
    polynomes_press = np.append(polyn_press, np.array([0, P0]))

    polyn_grav = np.polyfit(Radius[1:], Gravity[1:]/Radius[1:], 0)
    polynomes_gravity = np.append(polyn_grav, np.array([0]))
    
    maxlength = len(polynomes_temp)
    polynomes_all = np.zeros((maxlength, 3))
    polynomes_all[maxlength-len(polynomes_temp):maxlength,0] = polynomes_temp
    polynomes_all[maxlength-len(polynomes_press):maxlength,1] = polynomes_press
    polynomes_all[maxlength-len(polynomes_gravity):maxlength,2] = polynomes_gravity
    text = "Polynomial coefficients pour Temperature (K), Pressure (Pa), Gravity (m/s**2) as function of radius (m). Highest power first."
    np.savetxt(output_file1, polynomes_all, header=text, fmt='%.4e')
    print "Output file with polynomial coefficients: ", output_file1
    print "Number of coefficients (total): ", maxlength
    print "Number of free coefficients: ", Npoly+1
    print "==="
    
    r = np.linspace(Radius[0], Radius[-1], Npoints)
    temp = np.polyval(polyn_temp, r)*r**2+T0
    press = np.polyval(polyn_press, r)*r**2+P0
    
    ax[0,0].plot(r/1e3, temp, 'g')
    ax[1,0].plot(r/1e3, press, 'g')
    ax[0,1].plot(r/1e3, np.polyval(polyn_temp, r), 'g')
    ax[1,1].plot(r/1e3, np.polyval(polyn_press, r), 'g')

    

    ## polyn_temp = np.polyfit(Radius[1:], (Temperature[1:] -(T0))/Radius[1:]**2., 1)
    ## polyn_press = np.polyfit(Radius[1:], (Pressure[1:] -(P0))/Radius[1:]**2., 1)
    ## polynomes_temp = np.append(np.array([T0, 0]), polyn_temp)
    ## polynomes_press = np.append(np.array([P0, 0]), polyn_press)
    ## maxlength = len(polynomes_temp)
    ## polynomes_all = np.zeros((maxlength, 3))
    ## polynomes_all[0:len(polynomes_temp),0] = polynomes_temp
    ## polynomes_all[0:len(polynomes_press),1] = polynomes_press
    ## polynomes_all[0:len(polynomes_gravity),2] = polynomes_gravity

    ## text = "Polynomial coefficients pour Temperature (K), Pressure (Pa), Gravity (m/s**2) as function of radius (m). Highest power first."
    ## np.savetxt(output_file3, polynomes_all, header=text)
    ## print "Output file with polynomial coefficients - 2coeffs only : ", output_file1
    ## print "Number of coefficients: ", maxlength
    ## print "==="

    #new radius repartition (resampling with constant dr)
    r = np.linspace(Radius[0], Radius[-1], Npoints)
    temp = np.polyval(polyn_temp, r)*r**2+T0
    press = np.polyval(polyn_press, r)*r**2+P0
    grav = np.polyval(polyn_grav, r)*r
    
    ax[0,0].plot(r/1e3, temp, 'b')
    ax[1,0].plot(r/1e3, press, 'b')
    ax[2,0].plot(r/1e3, grav, 'b')
    ax[0,1].plot(r/1e3, np.polyval(polyn_temp, r), 'k')
    ax[1,1].plot(r/1e3, np.polyval(polyn_press, r), 'k')
    ax[2,1].plot(r/1e3, np.polyval(polyn_grav, r), 'k')
    ax[2,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    
    if resample == 1:
        data_resample = np.concatenate((np.array([r]).T, np.array([temp]).T, np.array([press]).T, np.array([grav]).T), axis=1)
        np.savetxt(output_file2, data_resample, header="radius (m), temperature (K)", fmt='%.4e')
        print "Output file with resample data ", output_file2
        print "Number of points in radius: ", Npoints
        print "==="


    print "Figures:" 
    print "red '+': data points, blue line: fit and resample"
    plt.show()
    
    return



if __name__ == '__main__':
    
    # name of the data file (default: earth.dat)
    if len(sys.argv)>1:
        Names = [sys.argv[1]]
    else: Names = ['earth', 'earth_0.5', 'mars', 'vesta', 'venus']

    Npoints = 200
    folder = "./dat/"
        
    for name in Names:
        homogeneous_sample_data(name, folder, Npoly=1)
