import glob
import os
import sys

import numpy as np


totals = {}
for flavor in "multimap_sim", "noise_sim", "obsmat":
    flavor_total = 0
    print(f"{flavor}")
    if flavor not in totals:
        totals[flavor] = {}
    fnames = sorted(glob.glob(f"{flavor}/times_*_*.txt"))
    for fname in fnames:
        parts = os.path.basename(fname).replace(".txt", "").split("_")
        #TELE = f"{parts[1]}_{parts[2]}"
        TELE = parts[2]
        band = parts[3]
        _, nodes, seconds = np.genfromtxt(fname, unpack=True)
        nhours = np.sum(nodes * seconds / 3600)
        if TELE not in totals[flavor]:
            totals[flavor][TELE] = {}
        if band in totals[flavor][TELE]:
            totals[flavor][TELE][band] += nhours
        else:
            totals[flavor][TELE][band] = nhours
        print(f"    {TELE} {band} {nhours:8.2f} node hours")

total_hours = 0
for flavor in totals:
    print(f"{flavor}")
    flavor_hours = 0
    for TELE in sorted(totals[flavor]):
        if TELE == "SPLAT":
            continue
        print(f"    {TELE}")
        tele_hours = 0
        for band in sorted(totals[flavor][TELE]):
            band_hours = totals[flavor][TELE][band]
            tele_hours += band_hours
            print(f"        {band}  {band_hours:10.2f}")
        flavor_hours += tele_hours
        print(f"        TOTAL {tele_hours:10.2f}")
    total_hours += flavor_hours
    print(f"    {flavor} TOTAL {flavor_hours:10.2f}")
print(f"TOTAL {total_hours:10.2f} total node hours")
