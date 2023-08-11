#!/usr/bin/env python

from glob import glob
import os
import sys

import numpy as np

telescopes = {
    "chlat" : "LAT0_CHLAT",
    "splat" : "LAT2_SPLAT",
    "spsat1" : "SAT1_SAT",
    "spsat2" : "SAT2_SAT",
    "spsat3" : "SAT3_SAT",
}

for tele, TELE in telescopes.items():
    if tele == "splat":
        freqs = [20, 30, 40, 90, 150, 220, 280]
    elif tele == "chlat":
        freqs = [30, 40, 90, 150, 220, 280]
    elif tele == "spsat1":
        freqs = [95, 155, 220, 280]
    elif tele == "spsat2":
        freqs = [85, 95, 145, 155, 220, 280]
    elif tele == "spsat3":
        freqs = [30, 40, 85, 145]
    else:
        raise RuntimeError(f"Unknown tele: {tele}")

    print(f"\ntelescope = {tele}, TELESCOPE = {TELE}")

    tele_type = tele[:5]

    #for suffix in "", "_lowcomplexity", "_highcomplexity":
    for suffix in "",:
        #print(f"\nsuffix = {suffix}\n")

        hits_tot = 0
        hits_tot_staged = 0
        nsplit_tot = 0

        for subset in "upto2mm", "over2mm":
            fnames_in = sorted(glob(f"../split_schedules_1_{subset}/{tele_type}/*.txt"))
            nsplit = len(fnames_in)
            if nsplit == 0:
                continue
            print(f"\nsubset = {subset}\n")
            print(f"Found {nsplit} split schedule files")
            obs = set()
            for fname in fnames_in:
                obs.add(os.path.basename(fname).replace(".txt",""))

            for freq in freqs:
                freq_obs_cleared = set()
                for fname in sorted(glob(f"cleared_logs{suffix}/{TELE}/f{freq:03}/*.log")):
                    freq_obs_cleared.add(os.path.basename(fname).replace(".log", ""))
                hits_cleared = 0
                for ob in freq_obs_cleared:
                    if ob in obs:
                        hits_cleared += 1
                hits = hits_cleared
                print(f"{freq:3} GHz is {hits:5} / {nsplit:5} = {hits / nsplit * 100:7.3f}% complete")
                nsplit_tot += nsplit
                hits_tot += hits

        print(f"\nTOTAL {hits_tot:5} / {nsplit_tot:5} = {hits_tot / nsplit_tot * 100:7.3f}% complete")
