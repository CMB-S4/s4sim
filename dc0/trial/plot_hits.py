import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

from toast.pixels_io import read_healpix, write_healpix

plt.rc("font", family="serif", size=24)
#plt.rc("text", usetex=True)

chlat = read_healpix("hitmap_outputs_chlat/mapmaker_hits.h5", nest=True, dtype=float)
splat = read_healpix("hitmap_outputs_splat/mapmaker_hits.h5", nest=True, dtype=float)
spsat = read_healpix("hitmap_outputs_spsat_2/mapmaker_hits.h5", nest=True, dtype=float)

spsat = hp.ud_grade(spsat, 256, order_in="NEST", order_out="NEST", power=-2)

for h in chlat, splat, spsat:
    h[h == 0] = hp.UNSEEN

nrow, ncol = 1, 3
fig = plt.figure(figsize=[9 * ncol, 6 * nrow])
cmap = "magma"
satmax = None # 100000

hp.mollview(spsat[0], nest=True, title="SPSAT", sub=[nrow, ncol, 1], cbar=False, cmap=cmap, max=satmax, xsize=1600)
hp.graticule(22.5)
hp.mollview(splat[0], nest=True, title="SPLAT", sub=[nrow, ncol, 2], cbar=False, cmap=cmap, xsize=1600)
hp.graticule(22.5)
hp.mollview(chlat[0], nest=True, title="CHLAT", sub=[nrow, ncol, 3], cbar=False, cmap=cmap, xsize=1600)
hp.graticule(22.5)

fig.savefig("hitmaps_row.pdf")

nrow, ncol = 3, 1
fig = plt.figure(figsize=[9 * ncol, 6 * nrow])
cmap = "magma"

hp.mollview(spsat[0], nest=True, title="SPSAT", sub=[nrow, ncol, 1], cbar=False, cmap=cmap, max=satmax, xsize=1600)
hp.graticule(22.5)
hp.mollview(splat[0], nest=True, title="SPLAT", sub=[nrow, ncol, 2], cbar=False, cmap=cmap, xsize=1600)
hp.graticule(22.5)
hp.mollview(chlat[0], nest=True, title="CHLAT", sub=[nrow, ncol, 3], cbar=False, cmap=cmap, xsize=1600)
hp.graticule(22.5)

fig.savefig("hitmaps_col.pdf")

plt.show()
