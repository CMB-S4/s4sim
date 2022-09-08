import os
import sys

import healpy as hp
import numpy as np

# Scalings as of 08/11/2022

# Which configurations go with which alternatives

alternatives = {
    "alternative_1" : [
        ("spsat", "spsat_baseline_deep"),
        ("spsat_aggressive", "spsat_aggressive_deep"),
    ],
    "alternative_2" : [
        ("spsat", "spsat_baseline_deep"),
        ("spsat_aggressive", "spsat_aggressive_deep"),
    ],
    "alternative_3" : [
        ("chsat_so", "chsat_baseline_so"),
        ("chsat_so_hwp", "chsat_baselinehwp_so"),
        ("chsat_so_hwp_aggressive", "chsat_aggressivehwp_so"),
        ("chsat_s4", "chsat_baseline_s4"),
        ("chsat_s4_hwp", "chsat_baselinehwp_s4"),
        ("chsat_s4_hwp_aggressive", "chsat_aggressivehwp_s4"),
    ],
}

# Number of optics tubes in each alternative

tube_counts = {
    "alternative_1" : {
        30 :    1,
        40 :    1,
        85 :    3,
        95 :    3,
        145 :   3,
        155 :   3,
        220 :   2,
        280 :   2,
    },
    "alternative_2" : {
        30 :    2,
        40 :    2,
        85 :    4,
        95 :    3,
        145 :   4,
        155 :   3,
        220 :   3,
        280 :   3,
    },
    "alternative_3" : {
        30 :    3,
        40 :    3,
        85 :    9,
        95 :    9,
        145 :   9,
        155 :   9,
        220 :   6,
        280 :   6,
    },
}

# Volume and NET scalings for a single optics tube

scalings = {
    "spsat" : {
        30 :    2.14,
        40 :    2.14,
        85 :   39.78,
        95 :   39.93,
        145 :  39.23,
        155 :  39.38,
        220 :  24.20,
        280 :  16.94,
    },
    "spsat_aggressive" : {
        30 :    (3.98, 2.138464718),
        40 :    (7.20, 1.181565162),
        85 :   (80.90, 1.521086954),
        95 :   (72.01, 1.561227106),
        145 :  (99.25, 1.222707841),
        155 :  (92.80, 1.194642016),
        220 :  (25.32, 0.954939689),
        280 :  (17.47, 0.9689532216),
    },
    "chsat_so" : {
        30 :    2.07,
        40 :    2.07,
        85 :   38.40,
        95 :   38.55,
        145 :  37.87,
        155 :  38.01,
        220 :  24.59,
        280 :  17.22,
    },
    "chsat_so_hwp" : {
        30 :    (1.58, 1.027741433),
        40 :    (1.58, 1.02767165),
        85 :   (29.29, 1.043875238),
        95 :   (28.84, 1.072304098),
        145 :  (28.31, 1.064990291),
        155 :  (28.69, 1.062919846),
        220 :  (19.15, 1.020298481),
        280 :  (13.25, 1.03233061),
    },
    "chsat_so_hwp_aggressive" : {
        30 :    (2.66, 2.109276308),
        40 :    (4.48, 1.251117391),
        85 :   (62.33, 1.897951595),
        95 :   (55.15, 1.869397564),
        145 :  (89.59, 1.302099193),
        155 :  (78.93, 1.28807286),
        220 :  (29.04, 1.030524203),
        280 :  (19.80, 1.058123266),
    },
    "chsat_s4" : {
        30 :    2.00,
        40 :    2.00,
        85 :   37.04,
        95 :   37.18,
        145 :  36.53,
        155 :  36.66,
        220 :  23.49,
        280 :  16.44,
    },
    "chsat_s4_hwp" : {
        30 :    (1.56, 1.027741433),
        40 :    (1.56, 1.02767165),
        85 :   (28.99, 1.043875238),
        95 :   (28.55, 1.072304098),
        145 :  (28.02, 1.064990291),
        155 :  (28.41, 1.062919846),
        220 :  (18.77, 1.020298481),
        280 :  (12.99, 1.03233061),
    },
    "chsat_s4_hwp_aggressive" : {
        30 :    (2.63, 2.109276308),
        40 :    (4.43, 1.251117391),
        85 :   (61.70, 1.897951595),
        95 :   (54.59, 1.869397564),
        145 :  (88.68, 1.302099193),
        155 :  (78.13, 1.28807286),
        220 :  (28.47, 1.030524203),
        280 :  (19.41, 1.058123266),
    },
}

for alternative, flavors in alternatives.items():
    n_optics_tube = tube_counts[alternative]
    for flavor_in, flavor_out in flavors:
        outdir = f"scaled_outputs/{alternative}/{flavor_out}"
        os.makedirs(outdir, exist_ok=True)

        for freq, scale in scalings[flavor_in].items():
            try:
                # The NET factor is included in `scale` but should not apply to hits
                scale, net_factor = scale
            except TypeError as e:
                net_factor = 1
            scale *= n_optics_tube[freq]
            base_flavor = flavor_in.replace("_hwp", "").replace("_aggressive", "")
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
