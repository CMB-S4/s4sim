import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


f_weights = {
    20 : 0.8,
    30 : 0.8,
    40 : 0.8,
    90 : 0.8,
    150 : 0.8,
    220 : 0.8,
    280 : 0.8,
}


ndets = {
    "phase1" : {
        20 : 216,
        30 : 768,
        40 : 768,
        90 : 2904 * 16,
        150 : 2904 * 16,
        220 : 1110 * 16,
        280 : 1110 * 16,
    },
    "phase2" : {
        20 : 216,
        30 : 768,
        40 : 768,
        90 : 2904 * 16,
        150 : 2904 * 16,
        220 : 1110 * 16,
        280 : 1110 * 16,
    },
}


def best3pc(h):
    hsorted = np.sort(h)
    limit = hsorted[int(0.40 * h.size)]
    return h >= limit


for freq in 20, 30, 40, 90, 150, 220, 280:
    band = f"f{freq:03}"

    f_weight = f_weights[freq]
    
    fname0 = f"../../phase1/sims/outputs/lat_wide/{band}/mapmaker_hits.fits"
    fname1 = f"outputs/lat_wide/{band}/season/mapmaker_hits.fits"
    fname2 = f"outputs/lat_wide/{band}/break/mapmaker_hits.fits"

    w1 = hp.read_map(fname0.replace("hits", "invcov"))
    w2season = hp.read_map(fname1.replace("hits", "invcov"))
    w2break = hp.read_map(fname2.replace("hits", "invcov"))
    w2 = w2season + w2break
    mask1 = best3pc(w1)
    mask2 = best3pc(w2)
    w1sum = np.sum(w1[mask1])
    w2sum = np.sum(w2[mask2])
                     
    h1 = hp.read_map(fname0)
    h2season = hp.read_map(fname1)
    h2break = hp.read_map(fname2)
    h2 = h2season + h2break

    h1sum = np.sum(h1[mask1])
    h2sum_season = np.sum(h2season[mask2])
    h2sum_break = np.sum(h2break[mask2])
    h2sum = np.sum(h2[mask2])

    var1 = np.sum(h1[mask1] / w1[mask1])
    var2 = np.sum(h2[mask2] / w2[mask2])

    ndet1 = ndets["phase1"][freq]
    ndet2 = ndets["phase2"][freq]

    print(
        f"{band} : "
        f" "
        f" survey weight = {w2sum / w1sum / f_weight:.2f}, "
        f" f_weight = {1 / f_weight:.2f}, "
        f" 1 / NET^2 = {var1 / var2:.2f}, "
        f" detector years = {h2sum / h1sum:.2f}, "
        f" detectors = {ndet2 / ndet1:.2f}, "
        f" tele years (season) = {(h2sum_season / ndet2) / (h1sum / ndet1):.2f}, x"
        f" tele years (off-season) = {h2sum / h2sum_season:.2f}"
    )
