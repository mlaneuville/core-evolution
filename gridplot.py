#!/usr/local/bin/python
# Time-stamp: <2016-03-01 15:09:53 marine>

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.mlab import griddata
import pandas as pd
import argparse
import sys

if __name__ == "__main__":

    basename = "out/"

    parser = argparse.ArgumentParser(description="Script to make a grid-plot from processed output.")
    parser.add_argument('-s', '--sub', type=str, help="subfolder name in which to search for convect.txt")
    parser.add_argument('-f', '--to_file', action='store_true', default=False, help="if set, will only write to file")

    args = parser.parse_args()

    if args.sub:
        basename += args.sub+"/"

    filename = basename+"/convect.txt"

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
    plt.savefig(basename+"convect-scatter.eps", format='eps', bbox_inches='tight')
    if not args.to_file:
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
    plt.savefig(basename+"convect-contourf.eps", format='eps', bbox_inches='tight')
    if not args.to_file:
        plt.show()
