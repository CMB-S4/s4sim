import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


if len(sys.argv) == 1:
    raise RuntimeError("Arguments: <covariance map> [<title>]")

fname_cov = sys.argv[1]
if len(sys.argv) > 2:
    title = sys.argv[2]
else:
    title = "fsky"


def fskies(hit):
    good = hit > 0
    npix = hit.size
    dhit = hit[good].astype(float)
    moment0 = np.sum(good) / npix
    moment1 = np.sum(dhit) / npix
    moment2 = np.sum(dhit ** 2) / npix
    moment4 = np.sum(dhit ** 4) / npix
    fraw = moment0
    fnoise = moment1 ** 2 / moment2
    fsignal = moment2 ** 2 / moment4
    return fraw, fnoise, fsignal


def invert_map(m):
    good = m != 0
    minv = np.zeros_like(m)
    minv[good] = 1 / m[good]
    return minv


def mask_map(m):
    masked = m.copy()
    masked[m == 0] = hp.UNSEEN
    return masked


def separate_hemispheres(m):
    """Construct Northern and Southern masks"""
    lon, lat = hp.pix2ang(hp.get_nside(m), np.arange(m.size), lonlat=True)
    rot = hp.Rotator(coord="cg")
    glon, glat = rot(lon, lat, lonlat=True)
    mask_south = glat < 0
    mask_north = glat > 0
    return mask_north, mask_south    


def plot_cov(cov, iplot, name, mask=None):
    invcov = invert_map(cov)
    rhit = invcov / np.amax(invcov)
    ftot = np.sum(rhit)

    if mask is not None:
        rhit *= mask
    fmask = np.sum(rhit)
    frac = fmask / ftot

    # Measure fsky

    fraw, fnoise, fsignal = fskies(rhit)
    
    # Convert to depth depth

    nside = hp.get_nside(cov)
    area = hp.nside2pixarea(nside)
    depth = np.sqrt(cov * area) * 1e6 * 180 / np.pi * 60  # uK.arcmin

    hp.mollview(
        mask_map(rhit),
        min=0,
        max=1,
        sub=[nrow, ncol, iplot],
        title=f"{name} fsky={fraw:.3f}/{fnoise:.3f}/{fsignal:.3f}, frac={frac:.3f}",
        unit="Relative",
        xsize=1600,
        cmap="magma",
    )


print(f"Loading {fname_cov}")
cov = hp.read_map(fname_cov, [3, 5])
cov = (cov[0] + cov[1]) / 2

print(f"Making hemisphere masks")
mask_north, mask_south = separate_hemispheres(cov)

print(f"Plotting")
nrow, ncol = 1, 3
fig = plt.figure(figsize=[ncol * 4, nrow * 3])
fig.suptitle(title)

plot_cov(cov, 1, "Full")
plot_cov(cov, 2, "South", mask_south)
plot_cov(cov, 3, "North", mask_north)

fname_plot = f"{title}.png"
fig.savefig(fname_plot)
print(f"Wrote {fname_plot}")

plt.show()
