#!/usr/bin/env python

from glob import glob
import os
import sys

import numpy as np


tele = "chlat"
TELE = "LAT0_CHLAT"

# logs/LAT0_CHLAT/f090/LAT0_CHLAT_split_schedule_0000_f090.log

#for suffix in "", "_lowcomplexity", "_highcomplexity":
for suffix in "",:
    print(f"\nsuffix = {suffix}\n")

    hits_tot = 0
    hits_tot_staged = 0
    nsplit_tot = 0

    for subset in "upto2mm", "over2mm":
        print(f"\nsubset = {subset}\n")
        fnames_in = sorted(glob(f"../split_schedules_1_{subset}/{tele}/*.txt"))
        #fnames_done = glob(f"split_schedules/{tele}/done/split_schedule_*.txt")
        #nsplit = len(fnames) + len(fnames_done)
        nsplit = len(fnames_in)
        print(f"Found {nsplit} split schedule files")
        obs = set()
        for fname in fnames_in:
            obs.add(os.path.basename(fname).replace(".txt",""))

        for freq in 30, 40, 90, 150, 220, 280:
            freq_obs_cleared = set()
            freq_obs_staged = set()
            for fname in sorted(glob(f"cleared_logs{suffix}/{TELE}/f{freq:03}/*.log")):
                freq_obs_cleared.add(os.path.basename(fname).replace(".log", ""))
            for fname in sorted(glob(f"staged_logs{suffix}/{TELE}/f{freq:03}/*.log")):
                freq_obs_staged.add(os.path.basename(fname).replace(".log", ""))
            hits_cleared = 0
            hits_staged = 0
            for ob in freq_obs_cleared:
                if ob in obs:
                    hits_cleared += 1
            for ob in freq_obs_staged:
                if ob in obs:
                    hits_staged += 1
            hits = hits_cleared + hits_staged
            print(f"{freq:3} GHz is {hits:5} / {nsplit:5} = {hits / nsplit * 100:7.3f}% complete", end="")
            print(f"  Staged {hits_staged:5} / {nsplit:5} = {hits_staged / nsplit * 100:7.3f}% complete")
            nsplit_tot += nsplit
            hits_tot += hits
            hits_tot_staged += hits_staged

    print(f"\nTOTAL {hits_tot:5} / {nsplit_tot:5} = {hits_tot / nsplit_tot * 100:7.3f}% complete")
    print(f"  Staged {hits_tot_staged:5} / {nsplit_tot:5} = {hits_tot_staged / nsplit_tot * 100:7.3f}% complete")
