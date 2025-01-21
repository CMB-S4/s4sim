import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


def invert_map(m):
    good = m != 0
    minv = np.zeros_like(m)
    minv[good] = 1 / m[good]
    return minv


def load_map(path):
    fname_season = f"{path}/season/mapmaker_cov.fits"
    fname_break = f"{path}/break/mapmaker_cov.fits"
    print(f"Loading {fname_season}")
    mseason = hp.read_map(fname_season)
    print(f"Loading {fname_break}")
    mbreak = hp.read_map(fname_break)
    return invert_map(invert_map(mseason) + invert_map(mbreak))


def plot_map(m, title, iplot, ref=None):
    # Transform pixel variance to depth
    mm = m.copy()
    nside = hp.get_nside(mm)
    area = hp.nside2pixarea(nside, degrees=True)
    mm = np.sqrt(mm * area) * 1e6 * 60  # uK.arcmin
    mm *= np.sqrt(2)  # I -> P
    good = m != 0
    mm[m == 0] = hp.UNSEEN
    mmsorted = np.sort(mm[good])
    limit = mmsorted[int(0.03 * mm.size)]
    best = np.logical_and(good, mm < limit)
    avg = np.mean(mm[best])
    title = title + f" {avg:.3f}"
    if ref is not None:
        title = title + f" (time : {(ref / avg)**2:.3f})"
    vmin = np.amin(mm[good])
    hp.mollview(
        mm,
        min=vmin,
        max=2 * vmin,
        cmap="magma",
        title=title,
        sub=[nrow, ncol, iplot],
        unit=r"$\mu$K.arcmin",
        format="%.3f",
    )
    return avg


band = "f155"

nrow, ncol = 2, 3
fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
#for iband, band in enumerate(bands):
m0 = load_map(f"outputs/sun90/{band}")
m90 = load_map(f"outputs/sun45_flag90/{band}")
m75 = load_map(f"outputs/sun45_flag75/{band}")
m60 = load_map(f"outputs/sun45_flag60/{band}")
m45 = load_map(f"outputs/sun45_flag45/{band}")
m1 = load_map(f"outputs/sun45/{band}")

avg = plot_map(m0, f"sun90", 1)
plot_map(m90, f"sun45/flag90", 2, avg)
plot_map(m75, f"sun45/flag75", 3, avg)
plot_map(m60, f"sun45/flag60", 4, avg)
plot_map(m45, f"sun45/flag45", 5, avg)
plot_map(m1, f"sun45", 6, avg)

fig.savefig("sat_depths_by_radius.png")
