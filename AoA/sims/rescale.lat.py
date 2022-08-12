import os
import sys

import healpy as hp
import numpy as np

# Scalings as of 08/05/2022

alternatives = {
    "alternative_1" : ["splat_for_spsat"],
    "alternative_2" : [
        "chlat_wide", "chlat_for_spsat", "chlat_for_spsat_w_wide", "chlat_ulf_for_spsat"
    ],
    "alternative_3" : [
        "chlat_for_chsat_so", "chlat_for_chsat_s4"
    ],
}

scalings = {
    "splat_for_spsat" : {
        20 :    85.98,
        30 :    98.27,
        40 :    98.27,
        90 :   405.95,
        150 :  378.52,
        220 :  374.30,
        280 :  258.71,
    },
    "chlat_wide" : {
        20 :     0.00,
        30 :   199.37,
        40 :   199.37,
        90 :   794.19,
        150 :  794.19,
        220 :  780.76,
        280 :  533.69,
    },
    "chlat_for_spsat" : {
        20 :     0.00,
        30 :   313.08,
        40 :   313.08,
        90 :  1247.19,
        150 : 1162.92,
        220 : 1140.95,
        280 :  788.59,
    },
    "chlat_for_spsat_w_wide" : {
        20 :     0.00,
        30 :   333.75,
        40 :   333.75,
        90 :  1329.55,
        150 : 1239.71,
        220 : 1216.29,
        280 :  840.67,
    },
    "chlat_for_chsat_so" : {
        20 :     0.00,
        30 :    98.48,
        40 :    98.48,
        90 :   392.31,
        150 :  392.31,
        220 :  385.68,
        280 :  263.63,
    },
    "chlat_for_chsat_s4" : {
        20 :     0.00,
        30 :    96.47,
        40 :    96.47,
        90 :   384.31,
        150 :  384.31,
        220 :  377.81,
        280 :  258.25,
    },
    "chlat_ulf_for_spsat" : {
        20 :  (104.31, 1.500660534),
        30 :   352.22,
        40 :   352.22,
        90 :  1247.19,
        150 : 1162.92,
        220 :  892.91,
        280 :  617.16,
    },
}

for alternative, flavors in alternatives.items():
    for flavor in flavors:
        outdir = f"scaled_outputs/{alternative}/{flavor}"
        os.makedirs(outdir, exist_ok=True)

        for freq, scale in scalings[flavor].items():
            try:
                # When net_factor is present, it will not be applied to the
                # hit maps.
                scale, net_factor = scale
            except TypeError as e:
                net_factor = 1
            if scale == 0:
                continue

            base_flavor = flavor.replace("_ulf", "")
            indir = f"outputs/lat/{base_flavor}/f{freq:03}"
            inmap = f"{indir}/mapmaker_hits.fits"
            outmap = f"{outdir}/hits_{freq:03}.fits"
            if freq == 20 and not os.path.isfile(inmap) and net_factor != 1:
                indir = f"outputs/lat/{base_flavor}/f030"
                inmap = f"{indir}/mapmaker_hits.fits"
                print(
                    f"WARNING: Using the 25GHz map as input for 20GHz output: "
                    f"{inmap} -> {outmap}"
                )
            print(f"Reading {inmap}")
            hits = hp.read_map(inmap, None)
            hits = (hits * scale * net_factor).astype(int)
            hp.write_map(outmap, hits, overwrite=True)
            print(f"Wrote {outmap}")

            inmap = f"{indir}/mapmaker_cov.fits"
            outmap = f"{outdir}/cov_{freq:03}.fits"
            print(f"Reading {inmap}")
            cov = hp.read_map(inmap, None)
            cov /= scale
            hp.write_map(outmap, cov, overwrite=True)
            print(f"Wrote {outmap}")
