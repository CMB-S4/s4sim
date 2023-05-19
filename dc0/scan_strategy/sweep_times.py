from glob import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np


if len(sys.argv) != 3:
    print("Usage: {sys.argv[0]} <schedule> <scan rate [deg/s]>")
    sys.exit()

fname = sys.argv[1]
scan_rate = float(sys.argv[2])

with open(fname, "r") as fin:
    fin.readline()
    header = fin.readline()
    try:
        site, telescope, lat, lon, alt = header.split()
    except:
        print(f"ERROR: failed to parse {fname}", flush=True)
        sys.exit()
    if np.abs(float(lat)) > 85:
        pole_site = True
    else:
        pole_site = False

# schedule = np.genfromtxt(fname, skip_header=3).T

arr = np.genfromtxt(
    fname,
    skip_header=3,
    usecols=[4, 5, 7, 8, 9, 10],
    dtype="f8,f8,S30,f8,f8,f8",
    names=["starts", "stops", "names", "az_mins", "az_maxs", "elevations"],
)
lengths = arr["stops"] - arr["starts"]
az_mins = arr["az_mins"]
az_maxs = arr["az_maxs"]
els = arr["elevations"]
throws = az_maxs - az_mins
throws[throws < 0] += 360
sweep_times = throws / scan_rate * np.cos(np.radians(els))
modulated_sweep_times = np.abs(np.cos(np.radians(az_mins)) - np.cos(np.radians(az_maxs))) / np.radians(scan_rate)

if lengths.size == 0:
    print(f"No scans found in {fname}")
    sys.exit()

if lengths.size == 1:
    print(f"Only one scan in {fname}")
    sys.exit()

good = throws != 0
if not np.all(good):
    print(f"WARNING: there are {np.sum(np.logical_not(good))} scans with zero throws")

print(
    f"Sweep time: {np.amin(sweep_times[good]):.1f} s "
    f"< {np.mean(sweep_times[good]):.1f} s "
    f"< {np.amax(sweep_times[good]):.1f} s"
)

print(
    f"Modulated sweep time: {np.amin(modulated_sweep_times[good]):.1f} s "
    f"< {np.mean(modulated_sweep_times[good]):.1f} s "
    f"< {np.amax(modulated_sweep_times[good]):.1f} s"
)
