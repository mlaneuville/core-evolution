#!/usr/bin/env python
# Time-stamp: "2015-12-07 16:34:59 marine"

# First argument : number of plots needed (1 if blank)
# Second argument : filenames are output-#.txt by default.
# If 'optional_run_name->output-XXX.txt' is used, add the optional_run_name as second argument.

import sys
import matplotlib.pyplot as plt
import numpy as np

basename = "output-"
end = ".txt"


if len(sys.argv)>1:
    N = int(sys.argv[1])
    print N, "files."
    if len(sys.argv)==3:
        basename = sys.argv[2]+"-"+basename
else: N = 0

fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    
for nombre in range(N+1):


    print "...", nombre
    rayon, temperature, conductivity = [], [], []


    with open(basename+str(nombre)+end,'r') as f:
        for i, l in enumerate(f):
            if len(l) > 4:
                L = l.split(" ")
                rayon = np.append(rayon, float(L[0]))
                temperature = np.append(temperature, float(L[2]))
                conductivity = np.append(conductivity, float(L[3]))
                if i==1:
                    Time = float(L[1])


    LABEL = "Time: "+"{:.0f}".format(Time)+" Ma"
    ax1.plot(rayon, temperature, label=LABEL)
    ax2.plot(rayon, conductivity, label=LABEL)

ax1.set_ylabel('Temperature')
ax2.set_ylabel('Conductivity')
ax2.set_xlabel('Radius')
fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
ax1.legend(loc=0)
plt.show()
