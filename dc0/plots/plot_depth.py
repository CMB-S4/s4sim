import os
import sys

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt


nrow, ncol = 3, 3
fig = plt.figure(figsize=[4 * ncol, 3 * nrow])
#tele = "chlat"
tele = "spsat"
freqs = {
    "chlat" : [25, 40, 90, 150, 230, 280],
    "splat" : [20, 25, 40, 90, 150, 230, 280],
    "spsat" : [25, 40, 85, 95, 145, 155, 230, 280],
}[tele]

for i, freq in enumerate(freqs):
    print(freq)
    # Cap the depth plot at MR 2.0
    if tele == "chlat":
        vmax = {
            25  : 21.8,
            40  : 12.4,
            90  :  2.0,
            150 :  2.0,
            230 :  6.9,
            280 : 16.7,
        }[freq]
    elif tele == "spsat":
        vmax = {
            25  :  3.5,
            40  :  4.5,
            85  :  0.88,
            95  :  0.78,
            145 :  1.2,
            155 :  1.3,
            230 :  3.5,
            280 :  6.0,
        }[freq]
        # translate Q/U requirement to I
        vmax /= np.sqrt(2)
    elif tele == "splat":
        vmax = {
            20  : 13.2,
            25  :  6.5,
            40  :  4.2,
            90  :  0.63,
            150 :  0.59,
            230 :  1.9,
            280 :  4.4,
        }[freq]
        # translate Q/U requirement to I
        vmax /= np.sqrt(2)
    else:
        raise RuntimeError(f"Unknown tele = {tele}")

    ivar = hp.read_map(
        f"/global/cfs/cdirs/cmbs4/dc/dc0/mission/{tele}/split01/{freq:03}/dc0_{tele}_t01.01_{freq:03}_mat02.fits"
    )

    nside = hp.get_nside(ivar)
    area = hp.nside2pixarea(nside, degrees=True)
    depth = np.sqrt(ivar * area) * 1e6 * 60
    depth[depth == 0] = hp.UNSEEN

    hp.mollview(
        depth,
        cmap="magma",
        max=vmax,
        title=f"{freq}GHz depth", unit="$\mu$K arcmin",
        sub=[nrow, ncol, i + 1],
    )
    hp.graticule(30)

plt.savefig(f"depth_{tele}.png")
plt.savefig(f"depth_{tele}.pdf")
