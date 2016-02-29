#!/usr/local/bin/python
# Time-stamp: <2016-02-25 21:22:13 marine>


import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":

    

    filename = "./out/convect.txt"

    dtype = "|S50" #"s40, s20, f8, f8, f8, f8, f8, f8, s20, f8" #  i8,f8,S5"# format
    data = np.genfromtxt(filename, dtype = dtype)

    a, b = data.shape
    print a, b

    convectstatus = data[:, 8]
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



    plt.scatter(TM.astype(float), diff.astype(float)*1e4*800, c = convection, s=100, cmap='prism')
    plt.xlabel("Mantle temperature (K)")
    plt.ylabel("Core conductivity (W/m/K")
    plt.show()
