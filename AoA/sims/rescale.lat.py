import os
import sys

import healpy as hp
import numpy as np

# Scalings as of 06/27/2022

scalings = {
    "splat_for_spsat" : {
        20 :    99.70,
        30 :    99.70,
        40 :    99.70,
        90 :   398.89,
        150 :  371.84,
        220 :  366.45,
        280 :  253.28,
    },
    "chlat_wide" : {
        30 :   199.37,
        40 :   199.37,
        90 :   797.46,
        150 :  797.46,
        220 :  787.49,
        280 :  538.29,
    },
    "chlat_for_spsat" : {
        30 :   204.20,
        40 :   204.20,
        90 :   816.79,
        150 :  761.60,
        220 :  750.56,
        280 :  518.77,
    },
    "chlat_for_spsat_w_wide" : {
        30 :   204.20,
        40 :   204.20,
        90 :   816.79,
        150 :  761.60,
        220 :  750.56,
        280 :  518.77,
    },
    "chlat_for_chsat_so" : {
        30 :   99.68,
        40 :   99.68,
        90 :  398.73,
        150 : 398.73,
        220 : 393.75,
        280 : 269.14,
    },
    "chlat_for_chsat_s4" : {
        30 :    97.65,
        40 :    97.65,
        90 :   390.59,
        150 :  390.59,
        220 :  385.71,
        280 :  263.65,
    },
    #"chlat_for_chsat_s4" : {
    #    30 :   105.54,
    #    40 :   105.54,
    #    90 :   422.17,
    #    150 :  393.65,
    #    220 :  387.94,
    #    280 :  268.14,
    #},
}

for flavor in scalings:
    outdir = f"scaled_outputs/lat/{flavor}"
    os.makedirs(outdir, exist_ok=True)

    for freq, scale in scalings[flavor].items():
        indir = f"outputs/lat/{flavor}/f{freq:03}"
        
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
