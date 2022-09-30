import os
import sys

import healpy as hp
import numpy as np

# Scalings as of 09/14/2022

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
        85 :    2.50,
        95 :    2.50,
        145 :   2.46,
        155 :   2.47,
        220 :   1.51,
        280 :   1.06,
    },
    "spsat_aggressive" : {
        30 :    (3.98, 2.138464718),
        40 :    (7.20, 1.181565162),
        85 :    (5.08, 1.521086954),
        95 :    (4.51, 1.561227106),
        145 :   (6.24, 1.222707841),
        155 :   (5.81, 1.194642016),
        220 :   (1.58, 0.954939689),
        280 :   (1.09, 0.9689532216),
    },
    "chsat_so" : {
        30 :    2.07,
        40 :    2.07,
        85 :    2.41,
        95 :    2.41,
        145 :   2.38,
        155 :   2.38,
        220 :   1.54,
        280 :   1.08,
    },
    "chsat_so_hwp" : {
        30 :    (1.58, 1.027741433),
        40 :    (1.58, 1.02767165),
        85 :    (1.83, 1.043875238),
        95 :    (1.81, 1.072304098),
        145 :   (1.78, 1.064990291),
        155 :   (1.80, 1.062919846),
        220 :   (1.20, 1.020298481),
        280 :   (0.83, 1.03233061),
    },
    "chsat_so_hwp_aggressive" : {
        30 :    (2.66, 2.109276308),
        40 :    (4.48, 1.251117391),
        85 :    (3.92, 1.897951595),
        95 :    (3.45, 1.869397564),
        145 :   (5.63, 1.302099193),
        155 :   (4.94, 1.28807286),
        220 :   (1.82, 1.030524203),
        280 :   (1.24, 1.058123266),
    },
    "chsat_s4" : {
        30 :    2.00,
        40 :    2.00,
        85 :    2.33,
        95 :    2.33,
        145 :   2.29,
        155 :   2.30,
        220 :   1.47,
        280 :   1.03,
    },
    "chsat_s4_hwp" : {
        30 :    (1.56, 1.027741433),
        40 :    (1.56, 1.02767165),
        85 :    (1.82, 1.043875238),
        95 :    (1.79, 1.072304098),
        145 :   (1.76, 1.064990291),
        155 :   (1.78, 1.062919846),
        220 :   (1.17, 1.020298481),
        280 :   (0.81, 1.03233061),
    },
    "chsat_s4_hwp_aggressive" : {
        30 :    (2.63, 2.109276308),
        40 :    (4.43, 1.251117391),
        85 :    (3.88, 1.897951595),
        95 :    (3.42, 1.869397564),
        145 :   (5.57, 1.302099193),
        155 :   (4.89, 1.28807286),
        220 :   (1.78, 1.030524203),
        280 :   (1.21, 1.058123266),
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
            indir = f"outputs_full/sat/{base_flavor}/f{freq:03}"

            inmap1 = f"{indir}/madam_hmap.fits"
            inmap2 = f"{indir}/mapmaker_hits.fits"
            outmap = f"{outdir}/hits_{freq:03}.fits"
            if os.path.isfile(inmap1):
                inmap = inmap1 
            else:
                inmap = inmap2
            print(f"Reading {inmap}")
            hits = hp.read_map(inmap, None)
            hits = (hits * scale * net_factor).astype(int)
            hp.write_map(outmap, hits, overwrite=True)
            print(f"Wrote {outmap}")

            inmap1 = f"{indir}/madam_wcov.fits"
            inmap2 = f"{indir}/mapmaker_cov.fits"
            if os.path.isfile(inmap1):
                inmap = inmap1
            else:
                inmap = inmap2
            outmap = f"{outdir}/cov_{freq:03}.fits"
            print(f"Reading {inmap}")
            cov = hp.read_map(inmap, None)
            cov[cov == hp.UNSEEN] = 0
            cov /= scale
            hp.write_map(outmap, cov, overwrite=True)
            print(f"Wrote {outmap}")
