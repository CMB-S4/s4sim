import os
import sys

import healpy as hp
import numpy as np

# Scalings as of 08/05/2022
# Still missing the "aggressive" scalings

scalings = {
    "spsat" : {
        30 :    4.29,
        40 :    4.29,
        85 :  238.69,
        95 :  239.60,
        145 : 235.37,
        155 : 236.27,
        220 :  96.78,
        280 :  67.75,
    },
    "spsat_aggressive" : {
        30 :    4.29,
        40 :    4.29,
        85 :  238.69,
        95 :  239.60,
        145 : 235.37,
        155 : 236.27,
        220 :  96.78,
        280 :  67.75,
    },
    "chsat_so" : {
        30 :    4.14,
        40 :    4.14,
        85 :  230.41,
        95 :  231.29,
        145 : 227.21,
        155 : 228.08,
        220 :  98.38,
        280 :  68.86,
    },
    "chsat_so_hwp" : {
        30 :    2.90,
        40 :    3.03,
        85 :  160.55,
        95 :  163.46,
        145 : 164.39,
        155 : 163.37,
        220 :  67.57,
        280 :  48.80,
    },
    "chsat_so_hwp_aggressive" : {
        30 :    2.90,
        40 :    3.03,
        85 :  160.55,
        95 :  163.46,
        145 : 164.39,
        155 : 163.37,
        220 :  67.57,
        280 :  48.80,
    },
    "chsat_s4" : {
        30 :    3.99,
        40 :    3.99,
        85 :  222.24,
        95 :  223.09,
        145 : 219.15,
        155 : 219.99,
        220 :  93.96,
        280 :  65.77,
    },
    "chsat_s4_hwp" : {
        30 :    2.87,
        40 :    3.00,
        85 :  158.93,
        95 :  161.82,
        145 : 162.73,
        155 : 161.72,
        220 :  66.23,
        280 :  47.83,
    },
    "chsat_s4_aggressive" : {
        30 :    2.87,
        40 :    3.00,
        85 :  158.93,
        95 :  161.82,
        145 : 162.73,
        155 : 161.72,
        220 :  66.23,
        280 :  47.83,
    },
}

for flavor in scalings:
    outdir = f"scaled_outputs/sat/{flavor}"
    os.makedirs(outdir, exist_ok=True)

    for freq, scale in scalings[flavor].items():
        base_flavor = flavor.replace("_hwp", "").replace("_aggressive", "")
        indir = f"outputs/sat/{base_flavor}/f{freq:03}"

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
