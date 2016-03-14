#!/usr/local/bin/python
# Time-stamp: <2016-03-01 15:09:53 marine>


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata


if __name__ == "__main__":

    

    filename = "./out/mantle_cond_10/convect.txt"

    dtype = "|S50" #"s40, s20, f8, f8, f8, f8, f8, f8, s20, f8" #  i8,f8,S5"# format
    data = np.genfromtxt(filename, dtype = dtype)

    a, b = data.shape
    print a, b

    convectstatus = data[:, 8]
    convecttime = data[:, 10]
    TM = data[:, 5]
    diff = data[:, 6]

    # print diff

    convection = []
    for a in convectstatus:
        if a == "convective":
            convection.append(1)
        elif a =="no-convection":
            convection.append(0)
        elif a =="transient":
            convection.append(2)
        else:
            print "problem"



    scatterplot = plt.scatter(TM.astype(float), diff.astype(float)*1e4*800, c = convecttime.astype(float), s=100, cmap='gnuplot')
    plt.xlabel("Mantle temperature (K)")
    plt.ylabel("Core conductivity (W/m/K)")
    plt.colorbar(scatterplot)
    plt.savefig("convect.eps", format='eps', bbox_inches='tight')
    plt.show()



    def grid(x, y, z, resX=100, resY=100):
        "Convert 3 column data to matplotlib grid"
        xi = np.linspace(min(x), max(x), resX)
        yi = np.linspace(min(y), max(y), resY)
        Z = griddata(x, y, z, xi, yi, interp="linear")
        X, Y = np.meshgrid(xi, yi)
        return X, Y, Z

    X, Y, Z = grid(TM.astype(float), diff.astype(float)*1e4*800, convecttime.astype(float))
    plt.contourf(X, Y, Z, 100)
    plt.show()
