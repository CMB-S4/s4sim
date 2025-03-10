import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


cov = hp.read_map("scaled_outputs/sun90max_f090_180years_cov.fits")
good = cov != 0
hit = np.zeros_like(cov)
hit[good] = 1 / cov[good]

nside = hp.get_nside(hit)
npix = hit.size
lon, lat = hp.pix2ang(nside, np.arange(npix), lonlat=True)

def fskies(hit):
    good = hit > 0
    dhit = hit[good].astype(float)
    moment0 = np.sum(good) / npix
    moment1 = np.sum(dhit) / npix
    moment2 = np.sum(dhit ** 2) / npix
    moment4 = np.sum(dhit ** 4) / npix
    fraw = moment0
    fnoise = moment1 ** 2 / moment2
    fsignal = moment2 ** 2 / moment4
    return fraw, fnoise, fsignal


nrow, ncol = 1, 3
fig = plt.figure(figsize=[ncol * 6, nrow * 4])

iplot = 0
for field in "all", "south", "north":
    iplot += 1
    if field == "all":
        hits = hit.copy()
    elif field == "south":
        hits = hit.copy()
        hits[np.logical_and(lon > 110, lon < 250)] = 0
    elif field == "north":
        hits = hit.copy()
        hits[np.logical_or(lon < 110, lon > 250)] = 0
    fraw, fnoise, fsignal = fskies(hits)
    hits /= np.amax(hits)
    ipix = np.argmax(hits)
    if True:
        clon = np.sum(lon * hits) / np.sum(hits)
        clat = np.sum(lat * hits) / np.sum(hits)
    else:
        clon = lon[ipix]
        clat = lat[ipix]
    hits[hits == 0] = hp.UNSEEN
    hp.mollview(
        hits,
        title=f"{field} : {fraw:.3f} {fnoise:.3f} {fsignal:.3f} : RA = {clon:.3f}, Dec = {clat:.3f}",
        sub=[nrow, ncol, iplot],
        cmap="magma",
        unit="Rhit",
        min=0,
        max=1,
    )
    hp.projplot(clon, clat, 'bx', lonlat=True)

fig.savefig("sat_fsky.png")
