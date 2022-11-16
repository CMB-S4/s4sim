#!/usr/bin/env python

from glob import glob
import os
import sys

import numpy as np


tele = "chlat"
TELE = "LAT0_CHLAT"

# logs/LAT0_CHLAT/f090/LAT0_CHLAT_split_schedule_0000_f090.log

hits_tot = 0
nsplit_tot = 0

for subset in "upto2mm", "over2mm":
    print(f"\nsubset = {subset}\n")
    fnames_in = sorted(glob(f"split_schedules_1_{subset}/{tele}/*.txt"))
    #fnames_done = glob(f"split_schedules/{tele}/done/split_schedule_*.txt")
    #nsplit = len(fnames) + len(fnames_done)
    nsplit = len(fnames_in)
    print(f"Found {nsplit} split schedule files")
    obs = set()
    for fname in fnames_in:
        obs.add(os.path.basename(fname).replace(".txt",""))

    for freq in 30, 40, 90, 150, 220, 280:
        freq_obs = set()
        fnames = sorted(glob(f"cleared_logs/{TELE}/f{freq:03}/*.log"))
        for fname in fnames:
            freq_obs.add(os.path.basename(fname).replace(".log", ""))
        hits = 0
        for ob in freq_obs:
            if ob in obs:
                hits += 1
        print(f"{freq:3} GHz is {hits:5} / {nsplit:5} = {hits / nsplit * 100:6.3f}% complete")
        nsplit_tot += nsplit
        hits_tot += hits

print(f"\nTOTAL {hits_tot:5} / {nsplit_tot:5} = {hits_tot / nsplit_tot * 100:6.3f}% complete")
nsplit_tot += nsplit
hits_tot += hits
