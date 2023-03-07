import os
import sys

import numpy as np


total = 0
for flavor in "cmb", "foreground", "noise":
    flavor_total = 0
    print(f"{flavor}")
    for band in "f030", "f040", "f090", "f150", "f220", "f280":
        fname = f"{flavor}_sim/times_{band}.txt"
        _, nodes, seconds = np.genfromtxt(fname, unpack=True)
        nhours = np.sum(nodes * seconds / 3600)
        flavor_total += nhours
        print(f"    {band} {nhours:8.2f} node hours")
    print(f"  {flavor_total:.2f} node hours")
    total += flavor_total
print(f"{total:.2f} node hours")
