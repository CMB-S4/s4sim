import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


sat_patches = {
    "shallow"      : [100, -41, -20, -67, 10],
    "intermediate" : [ 83, -45,  -3, -65,  1],
    "deep"         : [ 52, -52,  30, -60,  1],
}


# Plots targets and foregrounds

fig = plt.figure(figsize=[18, 6])
#plt.suptitle("fwhm = {} deg, nbin = {}".format(fwhm, nbin))

def plot_patches(patches, color="white", lw=2):
    for patch in patches:
        lonmin, latmax, lonmax, latmin, priority = patches[patch]
        lons, lats = [], []
        for lon in np.linspace(lonmin, lonmax, 100):
            lons.append(lon)
            lats.append(latmax)
        for lat in np.linspace(latmax, latmin, 100):
            lons.append(lonmax)
            lats.append(lat)
        for lon in np.linspace(lonmax, lonmin, 100):
            lons.append(lon)
            lats.append(latmin)
        for lat in np.linspace(latmin, latmax, 100):
            lons.append(lonmin)
            lats.append(lat)
        hp.projplot(lons, lats, '-', color="white", lw=4, threshold=1,
                    lonlat=True, coord='C')
        hp.projplot(lons, lats, '-', color="red", lw=2, threshold=1,
                    lonlat=True, coord='C')

def plot_p(fgmap, sub, title=None, nbin=10, coord="C", gr=15):
    """
    Plot the given map using only `nbin` distinct colors
    If `maps` contains several maps, the plotting color is
    chosen based on the most intense map: if any of the given
    """
    sorted_fgmap = np.sort(fgmap)
    lims = np.arange(nbin) / nbin
    npix = fgmap.size
    p = np.zeros(npix)
    for ilim, lim in enumerate(lims):
        val = sorted_fgmap[int(npix * lim)]
        p[fgmap > val] = lim
    hp.mollview(p, xsize=2400, sub=sub, title=title, coord=coord, cmap="inferno")
    hp.graticule(gr)
    return p

plot_p(hp.read_map("../fg_map_gal.fits"), sub=[1, 2, 1], coord="G")
plot_p(hp.read_map("../fg_map_equ.fits"), sub=[1, 2, 2], coord="C")
plot_patches(sat_patches, color="k", lw=3)

fig.savefig("fields_and_foregrounds.png")
