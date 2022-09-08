import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


year = 86400 * 365

for tele in "sat", "splat", "chlat":
    flavors = {
        "sat" : [
            "alternative_1/spsat_baseline_deep",
            "alternative_1/spsat_aggressive_deep", 
            "alternative_2/spsat_baseline_deep",
            "alternative_2/spsat_aggressive_deep",
            "alternative_3/chsat_baseline_s4",
            "alternative_3/chsat_baseline_so",
            "alternative_3/chsat_baselinehwp_s4",
            "alternative_3/chsat_baselinehwp_so",
            "alternative_3/chsat_aggressivehwp_s4",
            "alternative_3/chsat_aggressivehwp_so",
        ],
        "splat" : [
            "alternative_1/splat_tma_deep",
        ],
        "chlat" : [
            "alternative_1/chlat_cd_wide",
            "alternative_2/chlat_cd_hybrid",
            "alternative_2/chlat_cd_sp",
            "alternative_3/chlat_cd_s4",
            "alternative_3/chlat_cd_so",
            "alternative_3/chlat_cd_wide",
        ],
    }[tele]

    ndets = {
        "sat" : {
            "030" : 1 * 252,
            "040" : 1 * 252,
            "085" : 3 * 3084,
            "095" : 3 * 3552,
            "145" : 3 * 3084,
            "155" : 3 * 3552,
            "220" : 2 * 10008,
            "280" : 2 * 10008,
        },
        "splat" : {
            "020" : 216,
            "030" : 864,
            "040" : 864,
            "090" : 46656,
            "150" : 46656,
            "220" : 16884,
            "280" : 16884,
        },
        "chlat" : {
            "030" : 768,
            "040" : 768,
            "090" : 46656,
            "150" : 46656,
            "220" : 21574,
            "280" : 21574,
        }
    }

    if tele == "sat":
        fsample = 0.1
    else:
        fsample = 1.0

    for flavor in flavors:
        print(f"\n{flavor} :")
        print("  {:5}  {:9}  {:4}  {:7}  {:6}  {:6}  {:6}  {:10}  {:10}".format(
            "band", "telescope", "duty", "NET", "fsky", "fsky", "fsky", "I depth", "Q/U depth")
        )
        print("  {:5}  {:9}  {:4}  {:7}  {:6}  {:6}  {:6}  {:10}  {:10}".format(
            "[GHz]", "years", "[%]", "[uKrts]", "raw", "noise", "signal", "[uKarcmin]", "[uKarcmin]")
        )
        for band, ndet in ndets[tele].items():
            hit = hp.read_map(f"scaled_outputs/{flavor}/hits_{band}.fits")
            cov = hp.read_map(f"scaled_outputs/{flavor}/cov_{band}.fits")

            npix = hit.size
            good = hit > 0
            dhit = hit[good].astype(float)
            moment0 = np.sum(good) / npix
            moment1 = np.sum(dhit) / npix
            moment2 = np.sum(dhit ** 2) / npix
            moment4 = np.sum(dhit ** 4) / npix
            fraw = moment0
            fnoise = moment1 ** 2 / moment2
            fsignal = moment2 ** 2 / moment4

            """
            npix = hit.size
            fraw = np.sum(hit != 0) / npix
            fnoise = (np.sum(hit) / npix) ** 2 / (np.sum(hit ** 2) / npix)
            fsignal = (np.sum(hit ** 2) / npix) ** 2 / (np.sum(hit ** 4) / npix)
            """

            # Noise-weighted depth
            good = cov != 0
            covgood = cov[good].copy()
            #depth = covgood.size / np.sum(1 / covgood)
            depth = np.sum(1 / covgood) / np.sum(1 / covgood ** 2)
            #covgood = np.sort(covgood)
            #depth = np.mean(covgood[:covgood.size // 2])
            nside = hp.get_nside(cov)
            pixarea = hp.nside2pixarea(nside, degrees=True) * 60 ** 2
            depth = np.sqrt(depth * pixarea) * 1e6

            years = np.sum(hit) / ndet / fsample / year
            net = np.sqrt(np.median(cov[good] * hit[good]) / fsample) * 1e6
            print(
                f"  {band:>5}  {years:9.3f}  {years / 7 * 100:4.1f}"
                f"  {net:7.1f}"
                f"  {fraw:6.3f}  {fnoise:6.3f}  {fsignal:6.3f}"
                f"  {depth:10.3f}  {depth * np.sqrt(2):10.3f}",
            )
