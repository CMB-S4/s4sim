import os
import sys

import healpy as hp
import numpy as np

# Scalings as of 08/08/2022

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
        30 :    (7.96, 2.138464718),
        40 :   (14.40, 1.181565162),
        85 :  (485.41, 1.521086954),
        95 :  (432.06, 1.561227106),
        145 : (595.48, 1.222707841),
        155 : (556.80, 1.194642016),
        220 : (101.27, 0.954939689),
        280 :  (69.86, 0.9689532216),
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
        30 :    (3.15, 1.027741433),
        40 :    (3.15, 1.02767165),
        85 :  (175.74, 1.043875238),
        95 :  (173.07, 1.072304098),
        145 : (169.86, 1.064990291),
        155 : (172.17, 1.062919846),
        220 :  (76.60, 1.020298481),
        280 :  (53.00, 1.03233061),
    },
    "chsat_so_hwp_aggressive" : {
        30 :    (5.31, 2.109276308),
        40 :    (8.96, 1.251117391),
        85 :  (373.96, 1.897951595),
        95 :  (330.91, 1.869397564),
        145 : (537.52, 1.302099193),
        155 : (473.58, 1.28807286),
        220 : (116.18, 1.030524203),
        280 :  (79.20, 1.058123266),
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
        30 :    (3.12, 1.027741433),
        40 :    (3.12, 1.02767165),
        85 :  (173.96, 1.043875238),
        95 :  (171.32, 1.072304098),
        145 : (168.14, 1.064990291),
        155 : (170.43, 1.062919846),
        220 :  (75.08, 1.020298481),
        280 :  (51.95, 1.03233061),
    },
    "chsat_s4_hwp_aggressive" : {
        30 :    (5.26, 2.109276308),
        40 :    (8.87, 1.251117391),
        85 :  (370.19, 1.897951595),
        95 :  (327.57, 1.869397564),
        145 : (532.09, 1.302099193),
        155 : (468.80, 1.28807286),
        220 : (113.88, 1.030524203),
        280 :  (77.63, 1.058123266),
    },
}

for flavor in scalings:
    outdir = f"scaled_outputs/sat/{flavor}"
    os.makedirs(outdir, exist_ok=True)

    for freq, scale in scalings[flavor].items():
        try:
            # The NET factor is included in `scale` but should not apply to hits
            scale, net_factor = scale
        except TypeError as e:
            net_factor = 1
        base_flavor = flavor.replace("_hwp", "").replace("_aggressive", "")
        indir = f"outputs/sat/{base_flavor}/f{freq:03}"

        inmap = f"{indir}/mapmaker_hits.fits"
        outmap = f"{outdir}/hits_{freq:03}.fits"
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
