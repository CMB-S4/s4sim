import os
import sys

import numpy as np

for site in "chile", "pole":
    for tele in "sat", "lat":
        fname = f"{site}_{tele}/schedules/{site}_schedule_{tele}.txt"
        arr = np.genfromtxt(fname, skip_header=2).T
        times = arr[5] - arr[4]
        tsum = np.sum(times)
        print(f"{site} {tele} {tsum / 10}")
