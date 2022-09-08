import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


for alternative in os.listdir("scaled_outputs"):
    print(alternative)

    flavors = os.listdir(f"scaled_outputs/{alternative}")
    ncol = 3
    nrow = int(np.ceil(len(flavors) / ncol))
    fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
    fig.suptitle(f"{alternative} : Simulated depth at 150/155GHz")
    iplot = 0

    for flavor in sorted(flavors):
        print(alternative, flavor)
        if "LAT" in flavor.upper():
            band = "150"
        else:
            band = "155"
        iplot += 1
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
        vmin = np.amin(depth[good])
        sorted_depth = np.sort(depth[good])
        ngood = sorted_depth.size
        vmax = sorted_depth[int(0.75 * ngood)]
        hp.mollview(depth, min=vmin, max=vmax, title=f"{flavor}", sub=[nrow, ncol, iplot], unit="$\mu$K.arcmin")

    fig.savefig(f"depth_{alternative}.png")
