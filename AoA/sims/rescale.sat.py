import os
import sys

import healpy as hp
import numpy as np

# Scalings as of 06/27/2022 according to https://docs.google.com/spreadsheets/d/10_j-p4FmWuYgcpzr9Pvw6ZOQArhX2-OlDlvvkw9jd0w/edit#gid=1137632390

scalings = {
    "spsat" : {
        30 :    6.13,
        40 :    6.13,
        85 :  294.10,
        95 :  294.10,
        145 : 290.02,
        155 : 290.02,
        220 : 136.16,
        280 :  95.31,
    },
    "chsat_so" : {
        30 :    4.73,
        40 :    4.73,
        85 :  227.12,
        95 :  227.12,
        145 : 223.97,
        155 : 223.97,
        220 : 110.72,
        280 :  77.51,
    },
    "chsat_s4" : {
        30 :    4.56,
        40 :    4.56,
        85 :  219.07,
        95 :  219.07,
        145 : 216.02,
        155 : 216.02,
        220 : 105.75,
        280 :  74.02,
    },
}

for flavor in scalings:
    outdir = f"scaled_outputs/sat/{flavor}"
    os.makedirs(outdir, exist_ok=True)

    for freq, scale in scalings[flavor].items():
        indir = f"outputs/sat/{flavor}/f{freq:03}"
        
        inmap = f"{indir}/mapmaker_hits.fits"
        outmap = f"{outdir}/hits_{freq:03}.fits"
        print(f"Reading {inmap}")
        hits = hp.read_map(inmap, None)
        hits = (hits * scale).astype(int)
        hp.write_map(outmap, hits, overwrite=True)
        print(f"Wrote {outmap}")

        inmap = f"{indir}/mapmaker_cov.fits"
        outmap = f"{outdir}/cov_{freq:03}.fits"
        print(f"Reading {inmap}")
        cov = hp.read_map(inmap, None)
        cov /= scale
        hp.write_map(outmap, cov, overwrite=True)
        print(f"Wrote {outmap}")
