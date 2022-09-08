import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


for alternative in sorted(os.listdir("scaled_outputs")):
    print(alternative)

    flavors = os.listdir(f"scaled_outputs/{alternative}")
    ncol = 3
    nrow = int(np.ceil(len(flavors) / ncol))
    r = 0.75
    fig = plt.figure(figsize=[r * 4 * ncol, r * 4.0 * nrow])
    # fig.suptitle(f"{alternative} : Simulated depth at 150/155GHz")

    for flavor in sorted(flavors):
        print(alternative, flavor)
        if "LAT" in flavor.upper():
            band = "150"
        else:
            band = "155"
        iplot, title, vmin, vmax = {
            "alternative_1" : {
                "spsat_baseline_deep" :   (1, "SPSAT baseline",   0.5, 3.0),
                "spsat_aggressive_deep" : (2, "SPSAT aggressive", 0.5, 3.0),
                "splat_tma_deep" :        (3, "SPLAT delensing",  0.3, 2.4),
                "chlat_cd_wide" :         (4, "CHLAT wide",       1.2, 1.6),
            },
            "alternative_2" : {
                "spsat_baseline_deep"   : (1, "SPSAT baseline",   0.5, 3.0),
                "spsat_aggressive_deep" : (2, "SPSAT aggressive", 0.5, 3.0),
                "chlat_cd_sp" :           (3, "CHLAT delensing",  0.5, 3.0),
                "chlat_cd_wide" :         (4, "CHLAT wide",       1.2, 1.6),
                "chlat_cd_hybrid" :       (5, "CHLAT hybrid",     0.5, 3.0),
            },
            "alternative_3" : {
                "chsat_baseline_s4" :      (1, "CHSAT-S4 baseline",   0.5, 3.0),
                "chsat_baselinehwp_s4" :   (2, "CHSAT-S4 baseline+HWP",   0.5, 3.0),
                "chsat_aggressivehwp_s4" : (3, "CHSAT-S4 aggressive+HWP", 0.5, 3.0),
                "chsat_baseline_so" :      (4, "CHSAT-SO baseline",   0.5, 3.0),
                "chsat_baselinehwp_so" :   (5, "CHSAT-SO baseline+HWP",   0.5, 3.0),
                "chsat_aggressivehwp_so" : (6, "CHSAT-SO aggressive+HWP", 0.5, 3.0),
                "chlat_cd_s4" :            (7, "CHLAT-S4 delensing",  0.5, 3.0),
                "chlat_cd_so" :            (8, "CHLAT-SO delensing",  0.5, 3.0),
                "chlat_cd_wide" :          (9, "CHLAT wide",          1.2, 1.6),
            },
        }[alternative][flavor]
        vmin, vmax = 0.3, 2.4
        try:
            cov = hp.read_map(f"scaled_outputs/{alternative}/{flavor}/cov_{band}.fits")  # K^2
        except Exception as e:
            print(e)
            continue
        nside = hp.get_nside(cov)
        pix_area = hp.nside2pixarea(nside, degrees=True) * 60 ** 2  # arc min^2
        depth = np.sqrt(cov * pix_area) * 1e6  # uK arcmin
        good = depth != 0
        depth[depth == 0] = hp.UNSEEN
        #vmin = np.amin(depth[good])
        #sorted_depth = np.sort(depth[good])
        #ngood = sorted_depth.size
        #vmax = sorted_depth[int(0.75 * ngood)]
        hp.mollview(depth, min=vmin, max=vmax, title=title, sub=[nrow, ncol, iplot], unit="$\mu$K.arcmin", cmap="inferno")

    #fig.subplots_adjust(bottom=0.6)
    fig.savefig(f"depth_{alternative}.review.png")
