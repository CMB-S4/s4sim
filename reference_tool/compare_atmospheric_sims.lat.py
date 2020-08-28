import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


cmb = hp.read_map("out_v1/00000000/pole_cmb-unlensed_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits", nest=True)
atmo = hp.read_map("out_v1/00000000/pole_atmosphere_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits", nest=True)

good = atmo != hp.UNSEEN
atmo[good] *= 1e6

amp = 400
params = {
    #"rot" : [40, -52.5],
    "rot" : [40, -56.0],
    #"cmap" : "coolwarm",
    #"cmap" : "Blues_r",
    "cmap" : "bone",
    "xsize" : 1600,
    #"reso" : 1.50,
    "reso" : 0.70,
    "min" : -amp,
    "max" : amp,
    "nest" : True,
    "unit" : r"$\mu$K$_\mathrm{CMB}$",
    #"coord" : "C",
}

nrow, ncol = 2, 2
grparams = {
    "dpar" : 1,
    "dmer" : 15,
    #"coord" : "C",
    "local" : False,
}
plt.figure(figsize=[6 * ncol, 4 * ncol])

hp.gnomview(cmb, sub=[nrow, ncol, 1], **params, title="CMB")
hp.graticule(**grparams)

hp.gnomview(atmo, sub=[nrow, ncol, 2], **params, title="atmosphere, unscaled")
hp.graticule(**grparams)

# 2018 SPT mapmaking paper: 9087 hours of observing, 3491 independent maps, 150GHz observations
scale1 = (100 / 9087) ** .5
scale1 /= np.sqrt(1000 / 500) # We observe 1000 sq.deg, SPT saw 500.
scale2 = (10 / 3491) ** .5
atmo1 = atmo.copy()

atmo1[good] *= scale1
hp.gnomview(atmo1, sub=[nrow, ncol, 3], **params, title="atmosphere x {:.3f} (time)".format(scale1))
hp.graticule(**grparams)

atmo1[good] *= scale2 / scale1
hp.gnomview(atmo1, sub=[nrow, ncol, 4], **params, title="atmosphere x {:.3f} (visits)".format(scale2))
hp.graticule(**grparams)

plt.savefig("atmosphere_vs_spt.png")
