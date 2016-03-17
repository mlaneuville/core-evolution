#!/usr/local/bin/python
# Time-stamp: <2016-03-01 15:09:53 marine>

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.mlab import griddata
import pandas as pd
import sys

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "Enter the name of the subfolder where convect.txt is stored"
        print sys.argv[0]+" <folder>"
        sys.exit() 

    folder = sys.argv[1]
    filename = "./out/"+folder+"/convect.txt"

    names = ["code", "body", "TBL thickness", "TBL conductivity", "Kin viscosity", 
             "Mantle T", "Diff value", "Time max", "Status convection", "Onset convection",
             "Duration convection"]
    df = pd.read_table(filename, sep=' ', names=names, skiprows=1)
    dataset = df.sort(['Mantle T', 'Diff value'], ascending=[1, 1])
    print dataset.info()

    convectstatus = dataset["Status convection"]
    convecttime = dataset["Duration convection"]
    TM = dataset["Mantle T"]
    diff = dataset["Diff value"]

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

    # scatter plot
    scatterplot = plt.scatter(TM.astype(float), diff.astype(float)*1e4*800, c = convecttime.astype(float), s=100, cmap='jet')
    plt.xlabel("Mantle temperature (K)")
    plt.ylabel("Core conductivity (W/m/K)")
    plt.xlim(TM.values.min(), TM.values.max())
    plt.ylim(diff.values.min()*1e4*800, diff.values.max()*1e4*800)
    plt.colorbar(scatterplot)
    plt.savefig("convect1.eps", format='eps', bbox_inches='tight')
    plt.show()

    # contourf plot
    X, Y = np.meshgrid(TM.astype(float).unique(), diff.astype(float).unique()*1e4*800)
    Z = np.transpose(convecttime.astype(float).reshape(len(TM.unique()), len(diff.unique())))

    levels = [0, 250, 500, 750, 1000]

    CS1 = plt.contourf(X, Y, Z, levels, extend='min', cmap=cm.get_cmap("jet"))
    plt.colorbar(CS1, label="Convective era duration [Ma]")

    CS2 = plt.contour(X, Y, Z, levels, colors=('k',), linewidths=(3,))
    plt.clabel(CS2, fmt='%2.1f', colors='k', fontsize=14)

    plt.xlabel("Mantle temperature (K)")
    plt.ylabel("Core conductivity (W/m/K)")
    plt.savefig("convect2.eps", format='eps', bbox_inches='tight')
    plt.show()
