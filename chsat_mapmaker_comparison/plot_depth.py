import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

cov = hp.read_map("./outputs/noise/mapmaker_cov.fits")
nside = hp.get_nside(cov)
area = hp.nside2pixarea(nside, degrees=True)
depth = np.sqrt(cov) * np.sqrt(area) * 60 * 1e6
good = depth != 0
depth[depth == 0] = hp.UNSEEN
med = np.median(depth[good])
best = np.amin(depth[good])

hp.mollview(depth, max=med, cmap="magma", unit="$\mu$K arcmin", title="Depth")
hp.graticule(22.5)
plt.savefig("depth_full_sky.png")

imin = np.argmin(np.abs(depth))
lon, lat = hp.pix2ang(nside, imin, lonlat=True)

#hp.gnomview(depth, rot=[lon, lat], max=med, cmap="magma", unit="$\mu$K arcmin", title="Depth", xsize=800, reso=10)
rdepth = depth / best
rdepth[depth == hp.UNSEEN] = hp.UNSEEN
hp.gnomview(rdepth, rot=[int(lon), int(lat)], max=10, cmap="magma", unit="Max. depth", title="Relative depth", xsize=800, reso=10)
hp.graticule(5)
plt.savefig("depth_gnomview.png")

hp.orthview(rdepth, rot=[int(lon), int(lat)], max=10, cmap="magma", unit="Max. depth", title="Relative depth", xsize=800, half_sky=True)
hp.graticule(22.5)

npix = depth.size
print(np.sum(good) / npix)
print(np.sum(depth[good] < 3 * best) / npix)
plt.savefig("depth_orthview.png")
